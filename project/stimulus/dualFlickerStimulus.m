% @author: Radovan Vodila (radovan.vodila@ru.nl)
function flicker_protocol_two_images_hybrid
    %% ---- PARAMETERS ----
    devModeSkipSync = true;                % set false for real experiments
    flickerModeDefault = 'hybrid';         % freq, code, hybrid
    maxDisplaySec = 120;
    framesPerBit = 1;
    overlayAlphaDefault = 128;
    lb_lum = 50; hb_lum = 195;
    stimSize = 400;
    carrierHzs = [6, 12];
    imageFiles = {fullfile(pwd, 'project', 'stimulus', 'images', 'capybara.png'), ...
                  fullfile(pwd, 'project', 'stimulus', 'images', 'zebra2.png') ...
                  };
    ramp_len = 4;                          % frames for raised-cosine smoothing

    % Codes
    codefile = fullfile(pwd, 'project', 'stimulus', 'codes', 'mgold_61_6521.mat');
    S = load(codefile);
    code  = double(S.codes(1, :)); code2 = double(S.codes(2, :));
    code  = code(:)';            code2 = code2(:)';

    %% ---- DEFINE AREAS ----
    % INPUT (struct fields in 'args')
    %   attachStim    (int)    : 1..nStim; index of the base image this overlay attaches to.
    %   rel_x, rel_y  (0..1)   : fractional position inside the attached image rect
    %                            (0=left/top, 0.5=center, 1=right/bottom).
    %   w, h          (px)     : overlay rectangle size.
    %   alpha         (0..255) : overlay opacity used for alpha blending.
    %   lb, hb        (0..255) : luminance bounds; map01 → [lb, hb].
    %   flickerMode   (char)   : 'freq' | 'code' | 'hybrid'.
    %   carrierHz     (Hz)     : sine carrier frequency (used in 'freq'/'hybrid').
    %   code          (0/1 vec): binary code (used in 'code'/'hybrid'); expanded to frames.
    %   framesPerBit  (int≥1)  : frames per code bit.
    %   ramp_len      (int≥0)  : raised-cosine smoothing length at code transitions (frames).

    areas = [ ...
        makeArea(struct('attachStim',1,'rel_x',0.5,'rel_y',0.5,'w',stimSize,'h',stimSize, ...
                        'alpha',overlayAlphaDefault,'lb',lb_lum,'hb',hb_lum, ...
                        'flickerMode',flickerModeDefault,'carrierHz',carrierHzs(1), ...
                        'code',code,'framesPerBit',framesPerBit,'ramp_len',ramp_len)), ...
        makeArea(struct('attachStim',2,'rel_x',0.5,'rel_y',0.5,'w',stimSize,'h',stimSize, ...
                        'alpha',overlayAlphaDefault,'lb',lb_lum,'hb',hb_lum, ...
                        'flickerMode',flickerModeDefault,'carrierHz',carrierHzs(2), ...
                        'code',code2,'framesPerBit',framesPerBit,'ramp_len',ramp_len)), ...
    ];
    nAreas = numel(areas);

    %% ---- PSYCHTOOLBOX SETUP ----
    if devModeSkipSync
        Screen('Preference','SkipSyncTests', 1);  % dev: skip
    else
        Screen('Preference','SkipSyncTests', 0);  % real runs: enforce
    end 

    PsychDefaultSetup(2); KbName('UnifyKeyNames');
    screens = Screen('Screens'); screenNumber = max(screens);
    bgColor = [255 255 255];
    [win, winRect] = Screen('OpenWindow', screenNumber, bgColor);
    
    % Alpha blending mixes a new pixel (the source) with what’s already drawn (the destination) using an opacity α in [0,1]
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;
    ifi = Screen('GetFlipInterval', win); 
    displayFPS = 1/ifi;

    %% ---- LOAD & PLACE IMAGES ----
    nStim = numel(imageFiles);
    textures = zeros(1, nStim);
    dstRects = zeros(4, nStim);

    % positions: two images at 30% and 70% width, centered vertically
    rel_xs = linspace(0.3, 0.7, nStim);
    rel_y  = 0.5;

    winW = winRect(3); winH = winRect(4);
    for k = 1:nStim
        img = imread(imageFiles{k});
        if size(img,3)==1, img = repmat(img, [1 1 3]); end
        sz = size(img); minDim = min(sz(1:2));
        rowStart = floor((sz(1)-minDim)/2)+1;
        colStart = floor((sz(2)-minDim)/2)+1;
        imgSq = img(rowStart:rowStart+minDim-1, colStart:colStart+minDim-1, :);
        imgSq = imresize(imgSq, [stimSize stimSize]);
        imgSq = im2uint8(mat2gray(imgSq));  % keep identical behavior
        textures(k) = Screen('MakeTexture', win, imgSq);
        cx = winW * rel_xs(k); cy = winH * rel_y;
        dstRects(:,k) = CenterRectOnPointd([0 0 stimSize stimSize], cx, cy);
    end
    % Preload after textures exist
    Screen('PreloadTextures', win);

    %% ---- TIMEBASE ----
    totalFrames = max(round(maxDisplaySec / ifi), 1);
    t = (0:totalFrames-1) / displayFPS;

    %% ---- PRECOMPUTE MODULATION (BY AREA) ----
    [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t);

    %% ---- PRECOMPUTE OVERLAY RECTS (BY AREA) ----
    overlayRects = compute_overlay_rects_for_images(areas, dstRects);  % asserts inside

    %% ---- STIMULUS LOOP (optimized, VBL-locked) ----
    Priority(MaxPriority(win));

    % Prime pipeline: first flip may return immediately, second settles timing
    Screen('Flip', win);
    vbl = Screen('Flip', win);

    vbls       = zeros(1, totalFrames);
    targetVBLs = nan(1, totalFrames);     % purely diagnostic (theoretical)
    alphas     = double([areas.alpha]);    % 1 x nAreas, cached
    colors     = zeros(4, numel(areas));   % preallocate 4 x nAreas

    try
        for frameCount = 1:totalFrames
            % Draw base images
            Screen('FillRect', win, bgColor);               % cheap clear
            Screen('DrawTextures', win, textures, [], dstRects);

            % Draw overlays (batched) without re-allocating matrices
            lumRow = all_mod_lum(:, frameCount).';
            colors(1:3,:) = repmat(lumRow,3,1);
            colors(4,:)   = alphas;

            Screen('FillRect', win, colors, overlayRects);

            % Flip (block on next VBL; keeps pipeline busy)
            vblPrev = vbl;
            vbl = Screen('Flip', win);

            % Log actual and theoretical (for plots only)
            vbls(frameCount)       = vbl;
            targetVBLs(frameCount) = vblPrev + ifi;         % diagnostics, not used for scheduling

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
function area = makeArea(args)
% ARGUMENTS
% attachStim   : integer 1..nStim
%     Index of the base image this overlay attaches to (into textures/dstRects).
%
% rel_x, rel_y : scalar in [0,1]
%     Position INSIDE the attached image rect, fractional coords (0=left/top,
%     0.5=center, 1=right/bottom). Used as the overlay rect’s anchor point.
%
% w, h         : pixels
%     Overlay rectangle width/height drawn over the image (modulated patch).
%
% alpha        : 0..255
%     Overlay opacity used by alpha blending (0=transparent, 255=opaque).
%
% lb, hb       : 0..255
%     Lower/upper luminance bounds. The per-frame modulation is mapped into
%     this range: L = lb + (hb - lb) * map01, where map01 ∈ [0,1].
%
% flickerMode  : 'freq' | 'code' | 'hybrid'
%     'freq'   : sinusoidal carrier only (map01 = carrier in [0,1]).
%     'code'   : binary code only (map01 = code in [0,1]).
%     'hybrid' : carrier .* bipolar(code) in [-1,1], then map to [0,1] via
%                (x+1)/2 before lb/hb mapping (preserves full contrast).
%
% carrierHz    : Hz (≥0)
%     Sine frequency for 'freq' and 'hybrid' modes.
%
% code         : row/vector of 0/1
%     Binary sequence for 'code' and 'hybrid'. Expanded to frames by
%     framesPerBit and tiled to the full duration.
%
% framesPerBit : integer ≥1
%     Number of video frames per code bit. Bit rate = refreshHz / framesPerBit.
%     Increase to give more timing headroom on heavy scenes.
%
% ramp_len     : frames ≥0
%     Raised-cosine smoothing length applied at 0↔1 transitions of the code.
%     Reduces transients and missed deadlines at sharp edges.

defaults = struct('attachStim',1,'rel_x',0.5,'rel_y',0.5, ...
                  'w',200,'h',200,'alpha',128,'lb',60,'hb',200, ...
                  'flickerMode','hybrid','carrierHz',3,'code',[], ...
                  'framesPerBit',1,'ramp_len',2);
area = defaults;
fn = fieldnames(args);
for i=1:numel(fn), area.(fn{i}) = args.(fn{i}); end
end

function [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t)
% Compute per-area modulation & luminance sequences for all frames.
T = numel(t); nAreas = numel(areas);
mod_signals = zeros(nAreas, T);
all_mod_lum  = zeros(nAreas, T);
code_long_all = cell(1,nAreas);

for k = 1:nAreas
  A = areas(k);

  % Expand & smooth code if needed
  code_long = [];
  if ~isempty(A.code)
    code_expanded = repelem(A.code(:).', A.framesPerBit);
    nrep = ceil(T / numel(code_expanded));
    code_long = repmat(code_expanded, 1, nrep);
    code_long = code_long(1:T);
    pad_val = code_long(1);
    code_long_padded = [repmat(pad_val,1,A.ramp_len), code_long];
    code_long_smoothed = raised_cosine_smooth(code_long_padded, A.ramp_len);
    code_long = code_long_smoothed(A.ramp_len+1:end);
  end
  code_long_all{k} = code_long;

  carrier = 0.5 + 0.5 * sin(2*pi*A.carrierHz*t); % [0,1]
  switch lower(A.flickerMode)
    case 'freq'
      mod_signal = carrier;        map01 = mod_signal;       % [0,1]
    case 'code'
      assert(~isempty(code_long), 'Area %d is code-mode but code missing.', k);
      mod_signal = code_long;      map01 = mod_signal;       % [0,1]
    case 'hybrid'
      assert(~isempty(code_long), 'Area %d is hybrid-mode but code missing.', k);
      mod_signal = carrier .* (2*code_long - 1);             % [-1,1]
      map01 = (mod_signal + 1)/2;                            % -> [0,1]
    otherwise
      error('Unknown flickerMode: %s', A.flickerMode);
  end

  mod_signals(k,:) = mod_signal;
  all_mod_lum(k,:) = A.lb + (A.hb - A.lb) * map01;

  areas(k).code_long = code_long;
  areas(k).mod_signal = mod_signal;
end
end

function overlayRects = compute_overlay_rects_for_images(areas, dstRects)
% Rect per area, positioned relative to its attached base image rect.
n = numel(areas);
overlayRects = zeros(4,n);
nStim = size(dstRects,2);
for k = 1:n
  assert(areas(k).attachStim>=1 && areas(k).attachStim<=nStim, ...
      'Area %d attachStim=%d is out of range for %d images.', ...
      k, areas(k).attachStim, nStim);
  s = areas(k).attachStim;
  r = dstRects(:,s);  % [left top right bottom] of the base image
  w = r(3)-r(1); h = r(4)-r(2);
  cx = r(1) + areas(k).rel_x * w;
  cy = r(2) + areas(k).rel_y * h;
  overlayRects(:,k) = CenterRectOnPointd([0 0 areas(k).w areas(k).h], cx, cy);
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

fprintf('\n=== Frame Timing Diagnostics ===\n');
fprintf('Mean actual interval: %.5f s (%.2f Hz), SD: %.5f ms\n', ...
    mean(actual_intervals), 1/mean(actual_intervals), std(actual_intervals)*1000);
fprintf('Mean scheduled interval: %.5f s (%.2f Hz)\n', ...
    mean(scheduled_intervals), 1/mean(scheduled_intervals));
fprintf('Mean abs. timing error: %.5f ms (SD: %.5f ms)\n', ...
    mean(abs(timing_error))*1000, std(timing_error)*1000);

figure('Name','Image Flicker Diagnostics','NumberTitle','off');
tl = tiledlayout(5, max(2,nAreas), 'TileSpacing','compact');

% Row 1: code sequences
for k = 1:nAreas
    nexttile;
    if ~isempty(code_long_all{k})
        stairs(1:numel(code_long_all{k}), code_long_all{k}, 'LineWidth', 1.1);
        ylim([-0.2 1.2]); title(sprintf('Area %d: code',k));
    else
        plot(nan); ylim([0 1]); title(sprintf('Area %d: code (none)',k));
    end
    grid on; xlabel('Frame'); ylabel('Code');
end

% Row 2: luminance
for k = 1:nAreas
    nexttile;
    plot(t, all_mod_lum(k,:), 'LineWidth', 1.1);
    fm = flickerModeList(min(k,numel(flickerModeList)));
    chz = carrierHzs(min(k,numel(carrierHzs)));
    title(sprintf('Area %d: lum (%s, %.2f Hz)',k,fm, chz));
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


function code_smooth = raised_cosine_smooth(code_long, ramp_len)
% Raised-cosine smoothing for 0/1 sequences at transitions.
    code_smooth = code_long;
    N = numel(code_long);
    if ramp_len <= 0, return; end
    w = 0.5 * (1 - cos(linspace(0, pi, ramp_len))); % rising edge
    for i = 2:N
        if code_long(i) ~= code_long(i-1)
            j0 = i-ramp_len+1; if j0 < 1, continue; end
            if code_long(i) == 1
                code_smooth(j0:i) = w;       % rise 0->1
            else
                code_smooth(j0:i) = 1 - w;   % fall 1->0
            end
        end
    end
end
