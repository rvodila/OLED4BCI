% =========================================================================
%  flicker_protocol_two_images_hybrid
%  ------------------------------------------------------------------------
%  Author:  Radovan Vodila (radovan.vodila@ru.nl)
%
%  PURPOSE
%  -------
%  Present k flickering images with luminance modulation (freq / code / hybrid). 
%  Images are drawn with *feathered* (Gaussian) alpha at the edges so they fade smoothly into the background with no
%  hard corners. 
%  A third square (optosensor "diode") is drawn without feathering for clean photodiode recordings.
%
%  WHAT WE DO
%  ---------------------
%  1) Loads & centers images; square-crops and resizes to a fixed size.
%  2) Builds an RGBA texture for each image with a Gaussian edge alpha mask
%     (feather) — RGB = image, A = mask (0 at edges → 255 at center).
%  3) Builds a matching overlay texture with the same feathered alpha that
%     is tinted each frame to the desired luminance (per-stim modulation).
%  4) Runs a VBL-locked stimulus loop and flips at every refresh.
%  5) Prints and plots diagnostics (timing & modulation).
%
%  IMPORTANT RENDERING NOTES
%  -------------------------
%  - Blending is set to GL_SRC_ALPHA / GL_ONE_MINUS_SRC_ALPHA, so the alpha
%    channel of the image and overlay textures handles the feathered edges.
%  - The optosensor (stims(3)) is *not* feathered: it uses a full-alpha
%    square so the photodiode sees a sharp, high-contrast patch.
%
%  TODO / EXTENSIONS
%  -----------------
%  - Gamma linearization (make lb/hb photometrically linear)
%  - Optosensor integration with external trigger logging
%  - Add on-screen calibration helpers (photodiode probe)
% =========================================================================
function flicker_protocol_two_images_hybrid

    %% --------------------------------------------------------------------
    %  PARAMETERS, reference these in stimulus construction; 77 onwards
    % ---------------------------------------------------------------------
    devModeSkipSync     = true;          % set false for real experiments
    maxDisplaySec       = 60;             % total presentation time (s)
    
    flickerModeDefault  = 'hybrid';      % 'freq' | 'code' | 'hybrid' (default fallback)
    framesPerBit        = 2;             % code upsampling (frames/bit)
    overlayAlphaDefault = 128;           % default overlay alpha (used via stims)
    lb_lum              = 50;            % low luminance (0..255)
    hb_lum              = 195;           % high luminance (0..255)
    stimSize            = 400;           % image size after crop/resize (px)
    ramp_len            = 2;             % frames for raised-cosine smoothing of codes
    trialTaperFrames    = [];            % per-stim Hann taper at trial edges (frames)

    % ---- Codes (two distinct m-sequences for demonstration) --------------
    codefile = fullfile(pwd, 'project', 'stimulus', 'codes', 'mgold_61_6521.mat');
    S = load(codefile);
    code  = double(S.codes(1, :));
    code2 = double(S.codes(2, :));

    % ---- Feathering parameters (spatial Gaussian fade to background) -----
    %  featherCoreFrac : fraction of half-size that remains fully opaque.
    %                    Example: 0.6 → inner 60% (of radius) full alpha.
    %  featherSigmaPx  : Gaussian sigma in pixels controlling edge softness.
    featherCoreFrac = 0.6;
    featherSigmaPx  = 30;

    % ---- Stimulus definitions --------------------------------------------
    %  file        : path to image
    %  x,y         : screen center for the image
    %  size        : final size (px)
    %  alpha       : overlay alpha amplitude (0..255)
    %  lb / hb     : luminance range for modulation (0..255)
    %  flickerMode : 'freq' | 'code' | 'hybrid'
    %  carrierHz   : carrier frequency for freq/hybrid modes
    %  code        : bit sequence for code/hybrid modes
    %  framesPerBit: temporal upsampling of code bits
    %  ramp_len    : smoothing (frames) at code transitions
    %  trialTaperFrames : optional Hann taper length at trial begin/end
    stims(1) = orderfields(makeStim(struct( ...
        'file', fullfile(pwd,'project','stimulus','images','capybara.png'), ...
        'x', 640,  'y', 540, 'size', 400, ...
        'alpha', 128, 'lb', 50, 'hb', 195, ...
        'flickerMode', 'freq', 'carrierHz', 65, ...
        'code', code, 'framesPerBit', 2, 'ramp_len', 2, 'trialTaperFrames',60)));

    stims(2) = orderfields(makeStim(struct( ...
        'file', fullfile(pwd,'project','stimulus','images','zebra2.png'), ...
        'x', 1120, 'y', 540, 'size', 400, ...
        'alpha', 128, 'lb', 50, 'hb', 195, ...
        'flickerMode', 'freq', 'carrierHz', 60, ...
        'code', code2, 'framesPerBit', 2, 'ramp_len', 2,'trialTaperFrames',60)));

    % Optosensor (photodiode) square: kept *sharp* (no feather) by design.
    stims(3) = orderfields(makeStim(struct( ...
        'file', fullfile(pwd,'project','stimulus','images','white.png'), ...
        'x', 50, 'y', 50, 'size', 200, ...
        'alpha', 255, 'lb', 0, 'hb', 255, ...
        'flickerMode', 'freq', 'carrierHz', 120, ...
        'code', code2, 'framesPerBit', 2, 'ramp_len', 2)));
    diodeidx = 3; % index of the optosensor patch in stims

    %% --------------------------------------------------------------------
    %  PSYCHTOOLBOX SETUP
    % ---------------------------------------------------------------------
    if devModeSkipSync
        Screen('Preference','SkipSyncTests', 1);  % dev: skip sync tests
    else
        Screen('Preference','SkipSyncTests', 0);  % real runs: enforce sync
    end

    PsychDefaultSetup(2); KbName('UnifyKeyNames');
    screenNumber = max(Screen('Screens'));
    bgColor = [255 255 255];                         % white background
    [win, winRect] = Screen('OpenWindow', screenNumber, bgColor);

    % Alpha blending: source over destination with source alpha
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;

    ifi = Screen('GetFlipInterval', win);
    displayFPS = 1/ifi;                              % nominal refresh (Hz)

    %% --------------------------------------------------------------------
    %  LOAD, CROP, RESIZE, & PLACE IMAGES  (+ FEATHERED ALPHA)
    % ---------------------------------------------------------------------
    nStim    = numel(stims);
    textures = zeros(1,nStim);                       % RGBA image textures
    dstRects = zeros(4,nStim);                       % destination rects
    areas    = struct([]);                           % per-stim modulation params

    for k = 1:nStim
        assert(exist(stims(k).file,'file')==2, 'File not found: %s', stims(k).file);
        img = imread(stims(k).file);
        if size(img,3)==1, img = repmat(img,[1 1 3]); end

        % --- Center square-crop & resize to target size -------------------
        sz = size(img);
        minDim = min(sz(1:2));
        r0 = floor((sz(1)-minDim)/2)+1; c0 = floor((sz(2)-minDim)/2)+1;
        imgSq = img(r0:r0+minDim-1, c0:c0+minDim-1, :);
        imgSq = imresize(imgSq, [stims(k).size stims(k).size]);
        imgSq = im2uint8(mat2gray(imgSq));          % ensure 0..255

        % --- Feathered alpha: 1 at center, Gaussian falloff to 0 at edges -
        alphaMask = makeFeatherMask(stims(k).size, featherCoreFrac, featherSigmaPx);

        % --- Build RGBA texture: RGB=image, A=feathered alpha -------------
        rgba = zeros(stims(k).size, stims(k).size, 4, 'uint8');
        rgba(:,:,1:3) = imgSq;
        rgba(:,:,4)   = alphaMask;

        textures(k)   = Screen('MakeTexture', win, rgba);
        dstRects(:,k) = CenterRectOnPointd([0 0 stims(k).size stims(k).size], stims(k).x, stims(k).y);

        % --- Build matching "area" descriptor used for temporal modulation
        S = stims(k);
        areas(k).w = S.size;  areas(k).h = S.size;
        areas(k).alpha = S.alpha; areas(k).lb = S.lb; areas(k).hb = S.hb;
        areas(k).flickerMode = S.flickerMode; areas(k).carrierHz = S.carrierHz;
        areas(k).code = S.code; areas(k).framesPerBit = S.framesPerBit; areas(k).ramp_len = S.ramp_len;

        % Optional per-stim trial taper (Hann over begin/end of trial)
        if isfield(S,'trialTaperFrames') && ~isempty(S.trialTaperFrames)
            areas(k).trialTaperFrames = S.trialTaperFrames;
        else
            areas(k).trialTaperFrames = [];
        end

        % Optional area name (fallback to file base name)
        if isfield(S,'name') && ~isempty(S.name)
            areas(k).name = S.name;
        else
            [~, base] = fileparts(S.file);
            areas(k).name = base;
        end
    end



    %% --------------------------------------------------------------------
    %  BUILD OVERLAY TEXTURES (feathered alpha; diode kept sharp)
    % ---------------------------------------------------------------------
    % We draw a separate overlay per stim and tint its RGB each frame to the
    % desired luminance. The alpha (feather) lives in the texture itself.
    overlayTex  = zeros(1, nStim);
    for k = 1:nStim
        sz = stims(k).size;

        if k == diodeidx
            % Optosensor patch: NO feathering → full alpha (sharp square)
            a = uint8(255 * ones(sz, sz));
        else
            % Same feather as the image; also respect per-stim base alpha
            a = makeFeatherMask(sz, featherCoreFrac, featherSigmaPx);
            a = uint8(double(a) .* (double(stims(k).alpha)/255));
        end

        rgba = zeros(sz, sz, 4, 'uint8');
        rgba(:,:,1:3) = 255;    % white; will be per-frame tinted to target luminance
        rgba(:,:,4)   = a;      % feathered (or full) alpha
        overlayTex(k) = Screen('MakeTexture', win, rgba);
    end
    overlayRects = dstRects;     % overlays are 1:1 with image rects
    Screen('PreloadTextures', win);

    %% --------------------------------------------------------------------
    %  STIMULUS LOOP (VBL-locked)
    % ---------------------------------------------------------------------
    Priority(MaxPriority(win));

    % Precompute timing grid & diagnostics buffers
    waitframes  = 1;                                % update every refresh
    displayFPS  = 1/ifi;
    totalFrames = max(1, floor(maxDisplaySec * displayFPS / waitframes));
    t = (0:totalFrames-1) / displayFPS;             % seconds from start

    % Prime pipeline (first flip may return immediately; second stabilizes)
    Screen('Flip', win);
    vbl = Screen('Flip', win);
    targetVBLs = zeros(1, totalFrames);
    vbls       = zeros(1, totalFrames);

    % Precompute temporal modulation per area
    [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t);

    try
        for frameCount = 1:totalFrames
            % --- Draw base images ----------------------------------------
            Screen('FillRect', win, bgColor);                    % clear
            Screen('DrawTextures', win, textures, [], dstRects); % images (with feathered alpha)

            % --- Draw feathered overlays, tinted to per-frame luminance ---
            %     modColors: 4×n (RGB tint; A=255 so the texture's alpha masks the edges)
            lumRow    = all_mod_lum(:, frameCount).';                    % 1×nStim, in 0..255
            modColors = [repmat(lumRow,3,1); 255*ones(1,nStim)];         % RGB=tint, A=255
            Screen('DrawTextures', win, overlayTex, [], overlayRects, [], [], [], modColors);

            % --- Flip (locked to VBL) ------------------------------------
            vblPrev = vbl;
            vbl     = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi);

            vbls(frameCount)       = vbl;                                % actual flip time
            targetVBLs(frameCount) = vblPrev + waitframes * ifi;         % scheduled target

            % --- Lightweight exit check ---------------------------------
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

    %% --------------------------------------------------------------------
    %  DIAGNOSTICS (timing & modulation)
    % ---------------------------------------------------------------------
    flickerModes = string({areas.flickerMode});
    plot_modulation_diagnostics_img(mod_signals, all_mod_lum, t, code_long_all, ...
        flickerModes, vbls, targetVBLs, [areas.carrierHz], ifi, areas);
end

%% ===========================  HELPERS  ==================================

function [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t)
% PRECOMPUTE_AREA_MODULATIONS
% ---------------------------
% Build per-area temporal modulation signals for all frames.
% Returns:
%   areas(k).mod_signal  : map01 (0..1) sequence used for luminance mapping
%   mod_signals          : matrix [nAreas×T] of map01 sequences
%   all_mod_lum          : matrix [nAreas×T] of luminance values (0..255)
%   code_long_all        : expanded code (for diagnostics)
%
% Notes:
% - 'freq'    mode uses a sinusoidal carrier mapped around mid-gray.
% - 'code'    mode uses a (smoothed) binary code centered on 0.5.
% - 'hybrid'  mode multiplies bipolar code × bipolar carrier, then maps to 0..1.
T = numel(t); nAreas = numel(areas);
mod_signals   = zeros(nAreas, T);   % store map01 (0..1) for all modes
all_mod_lum   = zeros(nAreas, T);
code_long_all = cell(1,nAreas);

for k = 1:nAreas
    A = areas(k);

    % --- Trial-level Hann envelope (fade-in/out) --------------------------
    tf = 0;
    if isfield(A,'trialTaperFrames') && ~isempty(A.trialTaperFrames)
        tf = max(0, min(A.trialTaperFrames, floor(T/2)));  % clamp to ≤ T/2
    end
    env = ones(1,T);                        % default: no taper
    if tf > 0
        w = hann(2*tf)';                    % 0..1..0 across 2*tf points
        env(1:tf)           = w(1:tf);      % fade-in
        env(end-tf+1:end)   = w(tf+1:end);  % fade-out
    end

    % --- Expand & smooth code if needed ----------------------------------
    code_long = [];
    if ~isempty(A.code)
        code_expanded = repelem(A.code(:).', A.framesPerBit);
        nrep = ceil(T / numel(code_expanded));
        code_long = repmat(code_expanded, 1, nrep);
        code_long = code_long(1:T);

        % Safety warning if ramps too long for bit duration
        if A.ramp_len >= A.framesPerBit/2
            warning('ramp_len (%d) >= framesPerBit/2 (%g): ramps will overlap and flatten bits.', ...
                    A.ramp_len, A.framesPerBit/2);
        end

        code_long = raised_cosine_smooth(code_long, A.ramp_len); % smooth 0/1 edges
    end

    % --- Base carrier (map01: 0..1) --------------------------------------
    carrier01 = 0.5 + 0.5 * sin(2*pi*A.carrierHz*t);   % [0,1]

    % --- Mode selection with contrast-preserving taper -------------------
    switch lower(A.flickerMode)
        case 'freq'
            % Taper deviation around mid-gray to keep mean constant
            c = (carrier01 - 0.5) .* env;      % contrast in [-0.5,0.5]
            map01 = 0.5 + c;                   % back to [0,1]

        case 'code'
            assert(~isempty(code_long), 'Area %d is code-mode but code missing.', k);
            c = (code_long - 0.5) .* env;      % center 0/1 to [-0.5,0.5]
            map01 = 0.5 + c;                   % [0,1]

        case 'hybrid'
            assert(~isempty(code_long), 'Area %d is hybrid-mode but code missing.', k);
            % Bipolar code × bipolar carrier → taper → map to [0,1]
            mod_bipolar = (2*code_long - 1) .* (2*carrier01 - 1);  % [-1,1]
            mod_bipolar = mod_bipolar .* env;                      % taper
            map01 = (mod_bipolar + 1)/2;                           % [0,1]

        otherwise
            error('Unknown flickerMode: %s', A.flickerMode);
    end

    % Store
    mod_signals(k,:) = map01;                                  % [0,1]
    all_mod_lum(k,:) = A.lb + (A.hb - A.lb) * map01;           % luminance 0..255

    areas(k).code_long  = code_long;
    areas(k).mod_signal = map01;
    code_long_all{k}    = code_long;
end
end

function cleanup(win, textures)
% CLEANUP
% -------
% Close textures and the onscreen window; show cursor again.
    for i = 1:numel(textures), Screen('Close', textures(i)); end
    sca; ShowCursor;
end

function plot_modulation_diagnostics_img(mod_signals, all_mod_lum, t, code_long_all, flickerModeList, vbls, targetVBLs, carrierHzs, ifi, areas)
% PLOT_MODULATION_DIAGNOSTICS_IMG
% -------------------------------
% Visual diagnostics for:
%   - code sequences
%   - luminance time series
%   - autocorrelation per area
%   - timing error histogram & intervals
%   - flip timing error over time
nAreas = size(mod_signals,1);
actual_intervals = diff(vbls);
scheduled_intervals = diff(targetVBLs);
timing_error = vbls - targetVBLs;

% Collect names with fallback
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

% Row 2: luminance time series
for k = 1:nAreas
    nexttile;
    plot(t, all_mod_lum(k,:), 'LineWidth', 1.1);
    fm = flickerModeList(min(k,numel(flickerModeList)));
    chz = carrierHzs(min(k,numel(carrierHzs)));
    title(sprintf('%s — lum (%s, %.2f Hz)', names{k}, fm, chz));
    grid on; xlabel('Time (s)'); ylabel('Lum');
end

% % Row 3: autocorrelation per area
% for k = 1:nAreas
%     nexttile;
%     [acf, lags] = xcorr(mod_signals(k,:), 'coeff');
%     plot(lags, acf, 'LineWidth',1.1); ylim([0 1]);
%     title(sprintf('Area %d: autocorr',k));
%     grid on; xlabel('Lag (frames)'); ylabel('Norm. corr');
% end

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
% RAISED_COSINE_SMOOTH
% --------------------
% Smooth a 0/1 sequence using monotonic raised-cosine crossfades at each
% transition. ramp_len defines the half-window in samples on each side.
    y = x;
    N = numel(x);
    if ramp_len <= 0 || N < 3, return; end

    % Indices of NEW value (i+1 where diff~=0)
    trans = find(diff(x) ~= 0) + 1;

    for m = 1:numel(trans)
        i  = trans(m);                      % index where new value starts
        a0 = x(max(1, i-1));                % old level (0 or 1)
        a1 = x(i);                          % new level (0 or 1)

        i0 = max(1, i - ramp_len);          % left boundary
        i1 = min(N, i + ramp_len - 1);      % right boundary
        L  = i1 - i0 + 1;
        if L < 2, continue; end

        u  = linspace(0, 1, L);             % 0..1 across local window
        r  = 0.5 * (1 - cos(pi * u));       % monotonic raised-cosine 0→1
        y(i0:i1) = (1 - r) .* a0 + r .* a1; % crossfade old→new
    end
end

function s = makeStim(args)
% MAKESTIM
% --------
% Convenience constructor for the per-stimulus struct expected elsewhere.
% Fields:
%   file, x, y, size, alpha, lb, hb,
%   flickerMode, carrierHz, code, framesPerBit, ramp_len, trialTaperFrames
defaults = struct( ...
  'file','', 'x',[], 'y',[], 'size',400, ...
  'alpha',128,'lb',60,'hb',200, ...
  'flickerMode','hybrid','carrierHz',6,'code',[], ...
  'framesPerBit',1,'ramp_len',2, ...
  'trialTaperFrames',[]); % Hann taper frames (or [] for none)

s = defaults;
fn = fieldnames(args);
for i=1:numel(fn), s.(fn{i}) = args.(fn{i}); end
assert(~isempty(s.file) && ~isempty(s.x) && ~isempty(s.y), ...
  'Stim needs file, x, y.');
end

function a = makeFeatherMask(sz, coreFrac, sigmaPx)
% MAKEFEATHERMASK
% ---------------
% Create a uint8 alpha mask (0..255) with a fully-opaque circular core and
% Gaussian falloff to 0 at the edges.
% Inputs:
%   sz        : side length in pixels (mask is sz×sz)
%   coreFrac  : fraction of half-size that remains fully opaque (e.g., 0.6)
%   sigmaPx   : Gaussian sigma in pixels for the edge falloff
%
% The mask equals 255 inside radius r0 = coreFrac*(sz/2). Outside r0 the
% alpha decays as exp(-(r - r0)^2 / (2*sigma^2)).
    [X,Y] = meshgrid(1:sz, 1:sz);
    cx = (sz+1)/2; cy = (sz+1)/2;
    r  = sqrt((X-cx).^2 + (Y-cy).^2);
    r0 = coreFrac * (sz/2);                         % fully opaque inside r0
    g  = exp(-max(0, r - r0).^2 / (2*sigmaPx^2));  % Gaussian falloff outside
    a  = uint8(255 * min(max(g,0),1));
end
