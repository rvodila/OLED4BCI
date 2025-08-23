% @author: Radovan Vodila (radovan.vodila@ru.nl)
function flicker_protocol_two_images_hybrid
    %% ---- PARAMETERS ----
    devModeSkipSync     = true;          % set false for real experiments
    flickerModeDefault  = 'hybrid';      % 'freq' | 'code' | 'hybrid'
    maxDisplaySec       = 50;
    framesPerBit        = 2;
    overlayAlphaDefault = 128;
    lb_lum              = 50;
    hb_lum              = 195;
    stimSize            = 400;
    ramp_len            = 2;             % frames for raised-cosine smoothing
    trialTaperFrames    = [];   
    % Codes
    codefile = fullfile(pwd, 'project', 'stimulus', 'codes', 'mgold_61_6521.mat');
    S = load(codefile);
    code  = double(S.codes(1, :));
    code2 = double(S.codes(2, :));


    stims(1) = orderfields(makeStim(struct( ...
        'file', fullfile(pwd,'project','stimulus','images','capybara.png'), ...
        'x', 640,  'y', 540, 'size', 400, ...
        'alpha', 128, 'lb', 50, 'hb', 195, ...
        'flickerMode', 'freq', 'carrierHz', 6, ...
        'code', code, 'framesPerBit', 2, 'ramp_len', 2, 'trialTaperFrames',60)));

    stims(2) = orderfields(makeStim(struct( ...
        'file', fullfile(pwd,'project','stimulus','images','zebra2.png'), ...
        'x', 1120, 'y', 540, 'size', 400, ...
        'alpha', 128, 'lb', 50, 'hb', 195, ...
        'flickerMode', 'freq', 'carrierHz', 12, ...
        'code', code2, 'framesPerBit', 2, 'ramp_len', 2,'trialTaperFrames',60)));

    % OPTOSENSOR BOX (note: no stray quote after 255)
    stims(3) = orderfields(makeStim(struct( ...
        'file', fullfile(pwd,'project','stimulus','images','white.png'), ...
        'x', 50, 'y', 50, 'size', 200, ...
        'alpha', 250, 'lb', 0, 'hb', 255, ...
        'flickerMode', 'freq', 'carrierHz', 50, ...
        'code', code2, 'framesPerBit', 2, 'ramp_len', 2)));


    %% ---- PSYCHTOOLBOX SETUP ----
    if devModeSkipSync
        Screen('Preference','SkipSyncTests', 1);  % dev: skip
    else
        Screen('Preference','SkipSyncTests', 0);  % real runs: enforce
    end 

    PsychDefaultSetup(2); KbName('UnifyKeyNames');
    screenNumber = max(Screen('Screens'));
    bgColor = [255 255 255];
    [win, winRect] = Screen('OpenWindow', screenNumber, bgColor);
    
    % Alpha blending mixes a new pixel (the source) with what’s already drawn (the destination) using an opacity α in [0,1]
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;
    ifi = Screen('GetFlipInterval', win); 
    displayFPS = 1/ifi;

 %% ---- LOAD & PLACE IMAGES ----
nStim = numel(stims);
textures = zeros(1,nStim);
dstRects  = zeros(4,nStim);
areas     = struct([]);  % fresh

for k = 1:nStim
    assert(exist(stims(k).file,'file')==2, 'File not found: %s', stims(k).file);
    img = imread(stims(k).file);
    if size(img,3)==1, img = repmat(img,[1 1 3]); end

    % square-crop center + resize
    sz = size(img);
    minDim = min(sz(1:2));
    r0 = floor((sz(1)-minDim)/2)+1; c0 = floor((sz(2)-minDim)/2)+1;
    imgSq = img(r0:r0+minDim-1, c0:c0+minDim-1, :);
    imgSq = imresize(imgSq, [stims(k).size stims(k).size]);
    imgSq = im2uint8(mat2gray(imgSq));

    % make texture & rect
    textures(k) = Screen('MakeTexture', win, imgSq);
    dstRects(:,k) = CenterRectOnPointd([0 0 stims(k).size stims(k).size], stims(k).x, stims(k).y);

    % build matching area (once)
    S = stims(k);
    areas(k).w = S.size;  areas(k).h = S.size;
    areas(k).alpha = S.alpha; areas(k).lb = S.lb; areas(k).hb = S.hb;
    areas(k).flickerMode = S.flickerMode; areas(k).carrierHz = S.carrierHz;
    areas(k).code = S.code; areas(k).framesPerBit = S.framesPerBit; areas(k).ramp_len = S.ramp_len;

    % keep per-stim taper (don’t overwrite with a global [])
    if isfield(S,'trialTaperFrames') && ~isempty(S.trialTaperFrames)
        areas(k).trialTaperFrames = S.trialTaperFrames;
    else
        areas(k).trialTaperFrames = [];
    end

    % optional: area name
    if isfield(S,'name') && ~isempty(S.name)
        areas(k).name = S.name;
    else
        [~, base] = fileparts(S.file);
        areas(k).name = base;
    end
end

Screen('PreloadTextures', win);

%% overlays cover the images 1:1
overlayRects = dstRects;


    %% ---- STIMULUS LOOP (optimized, VBL-locked) ----
    Priority(MaxPriority(win));

    % Prime pipeline: first flip may return immediately, second settles timing
    waitframes  = 1;                    % 1 = update every refresh (use 2 to update every other)
    displayFPS  = 1/ifi;                % Hz
    % #flip iterations we will actually perform, honoring waitframes
    totalFrames = max(1, floor(maxDisplaySec * displayFPS / waitframes));
    t = (0:totalFrames-1) / displayFPS;
    
    Screen('Flip', win);          % prime
    vbl = Screen('Flip', win);    % prime
    targetVBLs = zeros(1, totalFrames);   % <= no NaNs
    vbls       = zeros(1, totalFrames);
    % preallocate
    alphas     = double([areas.alpha]);    % 1 x nAreas, cached
    colors     = zeros(4, numel(areas));   % preallocate 4 x nAreas
    % ----- PRECOMPUTE MODULATION (BY AREA) -----
    [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t);
    try
        for frameCount = 1:totalFrames
            % Draw base images
            Screen('FillRect', win, bgColor);               %  clear
            Screen('DrawTextures', win, textures, [], dstRects);

            % Draw overlays (batched) without re-allocating matrices
            lumRow = all_mod_lum(:, frameCount).';
            colors(1:3,:) = repmat(lumRow,3,1);
            colors(4,:)   = alphas;

            Screen('FillRect', win, colors, overlayRects);

            % Flip (VBL-locked)
            vblPrev = vbl;
            vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi);

            vbls(frameCount)       = vbl;
            targetVBLs(frameCount) = vblPrev + waitframes * ifi;  % current frame’s target

            % Lightweight exit check (optional: throttle to every few frames)
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE')), break; end
        end
    catch ME
        cleanup(win, textures);
        Priority(0);
        rethrow(ME);
    end
    cleanup(win, textures);
    Priority(0);

    %% ---- DIAGNOSTICS ----
    flickerModes = string({areas.flickerMode});
    plot_modulation_diagnostics_img(mod_signals, all_mod_lum, t, code_long_all, ...
        flickerModes, vbls, targetVBLs, [areas.carrierHz], ifi, areas);
end

%% ---------- HELPERS ----------


function [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t)
% Compute per-area modulation & luminance sequences for all frames.
T = numel(t); nAreas = numel(areas);
mod_signals   = zeros(nAreas, T);   % store map01 (0..1) for all modes
all_mod_lum   = zeros(nAreas, T);
code_long_all = cell(1,nAreas);

for k = 1:nAreas
    A = areas(k);

    % --- Build trial-level Hann envelope (0→1→0 over edges only) ---
    tf = 0;
    if isfield(A,'trialTaperFrames') && ~isempty(A.trialTaperFrames)
        tf = max(0, min(A.trialTaperFrames, floor(T/2)));  % clamp to ≤ T/2
    end
    env = ones(1,T);  % default: no taper
    if tf > 0
        w = hann(2*tf)';                 % 0..1..0 across 2*tf points
        env(1:tf)           = w(1:tf);   % fade-in
        env(end-tf+1:end)   = w(tf+1:end); % fade-out
    end

    % --- Expand & smooth code if needed ---
    code_long = [];
    if ~isempty(A.code)
        code_expanded = repelem(A.code(:).', A.framesPerBit);
        nrep = ceil(T / numel(code_expanded));
        code_long = repmat(code_expanded, 1, nrep);
        code_long = code_long(1:T);

        % SAFETY: ensure we don't obliterate bits
        if A.ramp_len >= A.framesPerBit/2
            warning('ramp_len (%d) >= framesPerBit/2 (%g): ramps will overlap and flatten bits.', ...
                    A.ramp_len, A.framesPerBit/2);
        end

        % Correct smoothing
        code_long = raised_cosine_smooth(code_long, A.ramp_len);
    end

    % --- Base carrier in [0,1] ---
    carrier01 = 0.5 + 0.5 * sin(2*pi*A.carrierHz*t);   % [0,1]

    % --- Mode selection with contrast-preserving tapering ---
    switch lower(A.flickerMode)
        case 'freq'
            % Taper deviation around mid-gray to keep mean constant
            c = (carrier01 - 0.5) .* env;      % contrast in [-0.5,0.5] tapered
            map01 = 0.5 + c;                   % back to [0,1]

        case 'code'
            assert(~isempty(code_long), 'Area %d is code-mode but code missing.', k);
            c = (code_long - 0.5) .* env;      % center 0/1 to [-0.5,0.5] then taper
            map01 = 0.5 + c;                   % [0,1]

        case 'hybrid'
            assert(~isempty(code_long), 'Area %d is hybrid-mode but code missing.', k);
            % Hybrid: bipolar product in [-1,1], taper contrast, then map to [0,1]
            mod_bipolar = (2*code_long - 1) .* (2*carrier01 - 1);  % [-1,1]
            mod_bipolar = mod_bipolar .* env;                      % taper
            map01 = (mod_bipolar + 1)/2;                           % [0,1]

        otherwise
            error('Unknown flickerMode: %s', A.flickerMode);
    end

    mod_signals(k,:) = map01;                                  % store [0,1]
    all_mod_lum(k,:) = A.lb + (A.hb - A.lb) * map01;           % luminance

    areas(k).code_long  = code_long;
    areas(k).mod_signal = map01;  
    code_long_all{k} = code_long;
end
end


function cleanup(win, textures)
    for i = 1:numel(textures), Screen('Close', textures(i)); end
    sca; ShowCursor;
end

function plot_modulation_diagnostics_img(mod_signals, all_mod_lum, t, code_long_all, flickerModeList, vbls, targetVBLs, carrierHzs, ifi, areas)
nAreas = size(mod_signals,1);
actual_intervals = diff(vbls); 
scheduled_intervals = diff(targetVBLs);
timing_error = vbls - targetVBLs;

% collect names with fallback
names = cell(1, nAreas);
for k = 1:nAreas
    if isfield(areas,'name') && numel(areas) >= k && ~isempty(areas(k).name)
        names{k} = areas(k).name;
    else
        names{k} = sprintf('Area %d', k);
    end
end

fprintf('\n=== Frame Timing Diagnostics ===\n');
fprintf('Mean actual interval: %.5f s (%.2f Hz), SD: %.5f ms\n', ...
    mean(actual_intervals), 1/mean(actual_intervals), std(actual_intervals)*1000);
fprintf('Mean scheduled interval: %.5f s (%.2f Hz)\n', ...
    mean(scheduled_intervals), 1/mean(scheduled_intervals));
fprintf('Mean abs. timing error: %.5f ms (SD: %.5f ms)\n', ...
    mean(abs(timing_error))*1000, std(timing_error)*1000);

figure('Name','Image Flicker Diagnostics','NumberTitle','off');
tl = tiledlayout(5, max(2,nAreas), 'TileSpacing','compact');


% Row 1: code sequences (binary -> stairs; smoothed -> plot)

for k = 1:nAreas
    nexttile;
    if ~isempty(code_long_all{k})
        stairs(1:numel(code_long_all{k}), code_long_all{k}, 'LineWidth', 1.1);
        ylim([-0.2 1.2]); title(sprintf('%s — code', names{k}));
    else
        plot(nan); ylim([0 1]); title(sprintf('%s — code (none)', names{k}));
    end
    grid on; xlabel('Frame'); ylabel('Code');
end


% Row 2: luminance
for k = 1:nAreas
    nexttile;
    plot(t, all_mod_lum(k,:), 'LineWidth', 1.1);
    fm = flickerModeList(min(k,numel(flickerModeList)));
    chz = carrierHzs(min(k,numel(carrierHzs)));
    title(sprintf('%s — lum (%s, %.2f Hz)', names{k}, fm, chz));
    grid on; xlabel('Time (s)'); ylabel('Lum');
end

% Row 3: autocorrelation per area
for k = 1:nAreas
    nexttile;
    [acf, lags] = xcorr(mod_signals(k,:), 'coeff');
    plot(lags, acf, 'LineWidth',1.1); ylim([0 1]);
    title(sprintf('Area %d: autocorr',k)); 
    grid on; xlabel('Lag (frames)'); ylabel('Norm. corr');
end

% Row 4: histogram + intervals
nexttile([1 max(1,ceil(nAreas/2))]);
histogram(timing_error*1000, 30); 
xlabel('Timing Error (ms)'); ylabel('Count'); 
title('Flip timing error histogram'); grid on;

nexttile([1 max(1,floor(nAreas/2))]);
h1 = plot(actual_intervals*1000,'-o'); hold on;
h2 = plot(scheduled_intervals*1000,'--');
h3 = yline(ifi*1000, 'k-');
ylabel('Frame Interval (ms)');
legend([h1 h2 h3], {'Actual','Scheduled',sprintf('IFI=%.2f ms',ifi*1000)}, 'Location','best');
grid on; title('Frame intervals');

% Row 5: time series of flip timing error
nexttile([1 max(1,nAreas)]); 
plot(timing_error*1000,'-o');
xlabel('Frame #'); ylabel('Timing error (ms)');
title('Flip timing error over time'); grid on;

nomHz = 1/ifi;
d = diff(vbls);
effHz = 1/mean(d);
jitter_ms = std(d)*1000;

title(tl, sprintf('Diagnostics | nAreas=%d | Nominal=%.2f Hz | Achieved=%.2f Hz | Jitter=%.2f ms', ...
                  nAreas, nomHz, effHz, jitter_ms));
end


function y = raised_cosine_smooth(x, ramp_len)
% Smooth 0/1 sequence x using monotonic raised-cosine crossfades at each transition.
% ramp_len = half-window in samples (on each side of the transition).
    y = x;
    N = numel(x);
    if ramp_len <= 0 || N < 3, return; end

    % indices of NEW value (i+1 where diff~=0)
    trans = find(diff(x) ~= 0) + 1;

    for m = 1:numel(trans)
        i  = trans(m);                % index where new value starts
        a0 = x(max(1, i-1));          % old level (0 or 1)
        a1 = x(i);                    % new level (0 or 1)

        i0 = max(1, i - ramp_len);    % left boundary
        i1 = min(N, i + ramp_len - 1);% right boundary (2*ramp_len samples total ideally)
        L  = i1 - i0 + 1;
        if L < 2, continue; end

        u  = linspace(0, 1, L);                 % 0..1 across the local window
        r  = 0.5 * (1 - cos(pi * u));           % monotonic raised-cosine 0->1
        y(i0:i1) = (1 - r) .* a0 + r .* a1;     % crossfade old->new
    end
end



%% 
function s = makeStim(args)
% IMAGE + COUPLED OVERLAY (overlay always centered, full image size)
% file      : image path
% x, y      : screen center (px) where the image goes
% size      : image size (px) after square-crop/resize
% alpha, lb, hb, flickerMode, carrierHz, code, framesPerBit, ramp_len

defaults = struct( ...
  'file','', 'x',[], 'y',[], 'size',400, ...
  'alpha',128,'lb',60,'hb',200, ...
  'flickerMode','hybrid','carrierHz',6,'code',[], ...
  'framesPerBit',1,'ramp_len',2, ...
  'trialTaperFrames',[]); % #frames for Hann taper (or [] for none) % fade length in frames (e.g., 96 by olaf)   

s = defaults;
fn = fieldnames(args);
for i=1:numel(fn), s.(fn{i}) = args.(fn{i}); end
assert(~isempty(s.file) && ~isempty(s.x) && ~isempty(s.y), ...
  'Stim needs file, x, y.');
end

