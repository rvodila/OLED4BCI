% @author: Radovan Vodila (radovan.vodila@ru.nl)
function flicker_protocol_video_hybrid

    %% ---- PARAMETERS ----
    devModeSkipSync = true;       % set false for real experiments
    flickerModeDefault = 'hybrid';% default per-area flickerMode
    overlayAlphaDefault = 128;
    lb_lum = 60; hb_lum = 200;
    framesPerBit = 1;
    carrierHzs = [3, 4];          % used in example areas below
    maxDisplayLen = 150;            % seconds (max playback duration)
    ramp_len = 4;                 % for raised cosine smoothing
    rectW = 300; rectH = 150;     % default area size
    movieFile = fullfile(pwd, 'project', 'stimulus', 'images', 'ape_walk.mp4'); % 25Hz, 17s, 950x540
    codefile = fullfile(pwd, 'project', 'stimulus', 'codes', 'mgold_61_6521.mat');
    audiofile = fullfile(pwd, 'project', 'stimulus', 'audio_files', 'chirpNoise_15chirps_0.68sDur_400Hz-4000Hz_17s.wav');

    % Load codes once
    S = load(codefile);
    code  = double(S.codes(1, :)); code2 = double(S.codes(2, :));
    code  = code(:)';              code2 = code2(:)';

    %% ---- DEFINE AREAS ----

    % MAKEAREA  Define a modulated overlay area for video presentation.
    %
    %   AREA = MAKEAREA(ARGS) returns a struct describing one overlay area.
    %   Fields in ARGS (all optional; defaults used if omitted):
    %
    %     rel_x, rel_y   Center position as fraction of video width/height (0–1).
    %     w, h           Width/height in pixels.
    %     alpha          Fill transparency (0=transparent, 255=opaque).
    %     lb, hb         Min/max luminance (0–255).
    %     flickerMode    'freq' (sine), 'code' (binary), or 'hybrid' (sine×code).
    %     carrierHz      Carrier frequency in Hz.
    %     code           0/1 vector for 'code'/'hybrid' modes; ignored for 'freq'.
    %     framesPerBit   Frames per code bit.
    %     ramp_len       Frames for raised-cosine smoothing at code transitions.
    %
    %   Used by modulation and geometry functions to render and modulate each area.

    areas = [ ...
        makeArea(struct('rel_x',0.4,'rel_y',0.3,'w',rectW,'h',rectH,'alpha',overlayAlphaDefault, ...
                        'lb',lb_lum,'hb',hb_lum,'flickerMode',flickerModeDefault,'carrierHz',carrierHzs(1), ...
                        'code',code,'framesPerBit',framesPerBit,'ramp_len',ramp_len)), ...
        makeArea(struct('rel_x',0.6,'rel_y',0.8,'w',rectW,'h',rectH,'alpha',overlayAlphaDefault, ...
                        'lb',lb_lum,'hb',hb_lum,'flickerMode',flickerModeDefault,'carrierHz',carrierHzs(2), ...
                        'code',code2,'framesPerBit',framesPerBit,'ramp_len',ramp_len)) ...
    ];
    nAreas = numel(areas);

    %% ---- AUDIO SETUP ----
    [audioData, audioFs] = audioread(audiofile);     % [samples x channels]
    nChannels = size(audioData,2);

    %% ---- PSYCHTOOLBOX SETUP ----
    if devModeSkipSync
      Screen('Preference','SkipSyncTests', 1);  % DEV: skip
    else
      Screen('Preference','SkipSyncTests', 0);  % EXPERIMENT: enforce
    end
    PsychDefaultSetup(2);
    KbName('UnifyKeyNames');


    screens = Screen('Screens');
    screenNumber = max(screens);
    bgColor = [128 128 128];
    [win, winRect] = Screen('OpenWindow', screenNumber, bgColor);
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;

    ifi = Screen('GetFlipInterval', win);
    displayFPS = 1/ifi;
    Priority(MaxPriority(win));

    %% ---- VIDEO SETUP ----
    videoReader = VideoReader(movieFile);
    videoFPS = videoReader.FrameRate;
    videoDuration = videoReader.Duration;
    nVidFrames = floor(videoFPS * videoDuration);

    % How many video frames fit in allotted time?
    nDisplayableFrames = min(round(videoFPS * maxDisplayLen), nVidFrames);
    actualDisplayLen   = nDisplayableFrames / videoFPS;

    % Timebases
    nDisplayFrames = round(actualDisplayLen * displayFPS); % # screen flips we will do
    totalFrames = nDisplayFrames;                          % single source of truth
    t = (0:totalFrames-1) / displayFPS;                    % seconds from first displayed frame

    % Build mapping from display frames to video frames (real-time playback)
    videoFrameTimes = (0:nDisplayableFrames-1) / videoFPS; % times of each video frame
    displayTimes    = (0:totalFrames-1) / displayFPS;
    videoFrameIdx   = interp1(videoFrameTimes, 1:nDisplayableFrames, displayTimes, 'previous', 1);
    videoFrameIdx   = min(max(round(videoFrameIdx),1), nVidFrames); % safety clamp
    useIdx          = unique(videoFrameIdx);               % frames we will actually display

    %% ---- PRECOMPUTE MODULATIONS PER AREA ----
    [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t);

    %% ---- LOAD ONLY NEEDED VIDEO FRAMES TO TEXTURES ----
    % We read frames sequentially and only MakeTexture for those in useIdx.
    % This saves GPU memory and a lot of MakeTexture calls.
    useMask = false(1, nVidFrames);
    useMask(useIdx) = true;

    videoTextures = containers.Map('KeyType','double','ValueType','double');
    for f = 1:nVidFrames
        frame = readFrame(videoReader);  % sequential read from disk
        if useMask(f)
            if size(frame,3)==1
                frame = repmat(frame, [1 1 3]);
            end
            videoTextures(f) = Screen('MakeTexture', win, im2uint8(frame));
        end
    end

    % Optional RAM printouts (Windows only)
    if ispc
        info = whos;
        totalBytes = sum([info.bytes]);
        fprintf('[Workspace] Total: %.2f MB (%.2f GB)\n', totalBytes/2^20, totalBytes/2^30);
    end

    %% ---- PRECOMPUTE GEOMETRY ----
    winW = winRect(3); winH = winRect(4);
    vidRectW = winW / 2;
    vidRectH = winH / 2;
    dstRect = CenterRectOnPointd([0 0 vidRectW vidRectH], winW/2, winH/2);

    overlayRects = compute_overlay_rects(areas, dstRect);

    %% ---- AUDIO DEVICE OPEN & FILL ----
    InitializePsychSound(1); % low-latency mode
    pahandle = PsychPortAudio('Open', [], 1, 1, audioFs, nChannels);
    PsychPortAudio('FillBuffer', pahandle, audioData'); % buffer must be [channels x samples]

    %% ---- VIDEO STIMULUS LOOP ----
    vbls       = zeros(1, totalFrames);
    targetVBLs = nan(1, totalFrames);   % diagnostics only

    try
        % Prime to real VBL cadence
        Screen('Flip', win);                  % first may return immediately
        vbl = Screen('Flip', win);            % aligned to VBL

        % Start audio at the NEXT VBL to phase-lock A/V
        PsychPortAudio('Start', pahandle, 1, vbl + ifi, 1);

        % Prealloc overlay colors
        alphas = double([areas.alpha]);       % 1 x nAreas
        colors = zeros(4, nAreas);            % 4 x nAreas

        for frameCount = 1:totalFrames
            % Clear background (if video doesn't cover full screen)
            Screen('FillRect', win, bgColor);

            % Draw video frame
            vidIdx = videoFrameIdx(frameCount);
            tex    = videoTextures(vidIdx);
            Screen('DrawTexture', win, tex, [], dstRect);

            % Batch overlays
            lumRow = all_mod_lum(:, frameCount).';
            colors(1:3,:) = repmat(lumRow,3,1);
            colors(4,:)   = alphas;
            Screen('FillRect', win, colors, overlayRects);

            % (Optional) signal end-of-draws
            % Screen('DrawingFinished', win, 0);

            % Flip at next VBL (unscheduled)
            vblPrev = vbl;
            vbl     = Screen('Flip', win);

            % Log actual & theoretical (for plots only)
            vbls(frameCount)       = vbl;
            targetVBLs(frameCount) = vblPrev + ifi;

            % ESC to abort (optional: throttle)
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE')), break; end
        end

    catch ME
        cleanup_all();
        rethrow(ME);
    end

    % Normal cleanup
    cleanup_all();


    %% ---- DIAGNOSTICS ----
    plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, ...
        flickerModeDefault, vbls, targetVBLs, [areas.carrierHz], ifi, areas);

    %% ---- NESTED CLEANUP TO CLOSE EVERYTHING ONCE ----
    function cleanup_all()
        try PsychPortAudio('Stop', pahandle); end %#ok<TRYNC>
        try PsychPortAudio('Close', pahandle); end %#ok<TRYNC>
        try closeTexturesMap(videoTextures); end %#ok<TRYNC>
        try cleanUpVideo(win); end %#ok<TRYNC>
        try Priority(0); end %#ok<TRYNC>
        ShowCursor;
    end
end

%% ---------- HELPERS ----------

function area = makeArea(args)
% Create a single overlay area with defaults.
defaults = struct( ...
  'rel_x',0.5, 'rel_y',0.5, ...
  'w',300, 'h',150, 'alpha',128, ...
  'lb',60, 'hb',200, ...
  'flickerMode','hybrid', ...     % 'freq'|'code'|'hybrid'
  'carrierHz',3, ...
  'code',[], ...                  % row vector of 0/1
  'framesPerBit',1, ...
  'ramp_len',2);
area = defaults;
fn = fieldnames(args);
for i=1:numel(fn)
  area.(fn{i}) = args.(fn{i});
end
end

function [areas, mod_signals, all_mod_lum, code_long_all] = precompute_area_modulations(areas, t)
% Precompute per-area modulation signals & luminance.
T = numel(t);
nAreas = numel(areas);
mod_signals = zeros(nAreas, T);
all_mod_lum  = zeros(nAreas, T);
code_long_all = cell(1,nAreas);

for k = 1:nAreas
  A = areas(k);

  % expand & smooth code if used
  code_long = [];
  if ~isempty(A.code)
    code_expanded = repelem(A.code(:).', A.framesPerBit);
    nrep = ceil(T / numel(code_expanded));
    code_long = repmat(code_expanded, 1, nrep);
    code_long = code_long(1:T);
    % smooth transitions with raised cosine
    pad_val = code_long(1);
    code_long_padded = [repmat(pad_val,1,A.ramp_len), code_long];
    code_long_smoothed = raised_cosine_smooth(code_long_padded, A.ramp_len);
    code_long = code_long_smoothed(A.ramp_len+1:end);
  end
  code_long_all{k} = code_long;

  % carrier always in [0,1]
  carrier = 0.5 + 0.5 * sin(2*pi*A.carrierHz*t);

  switch lower(A.flickerMode)
    case 'freq'
      mod_signal = carrier;          % [0,1]
      map01 = mod_signal;
    case 'code'
      if isempty(code_long), error('Area %d uses code mode but no code provided.', k); end
      mod_signal = code_long;        % [0,1]
      map01 = mod_signal;
    case 'hybrid'
      if isempty(code_long), error('Area %d uses hybrid mode but no code provided.', k); end
      mod_signal = carrier .* (2*code_long - 1); % [-1,1]
      map01 = (mod_signal + 1)/2;                % -> [0,1]
    otherwise
      error('Unknown flickerMode: %s', A.flickerMode);
  end

  mod_signals(k,:) = mod_signal;
  all_mod_lum(k,:) = A.lb + (A.hb - A.lb) * map01;

  % store for potential debugging
  areas(k).code_long = code_long;
  areas(k).mod_signal = mod_signal;
end
end

function overlayRects = compute_overlay_rects(areas, dstRect)
% Compute rects once from relative coordinates.
n = numel(areas);
overlayRects = zeros(4,n);
w = dstRect(3)-dstRect(1); h = dstRect(4)-dstRect(2);
for k = 1:n
  x = dstRect(1) + areas(k).rel_x * w;
  y = dstRect(2) + areas(k).rel_y * h;
  overlayRects(:,k) = CenterRectOnPointd([0 0 areas(k).w areas(k).h], x, y);
end
end

function closeTexturesMap(m)
% Close all textures stored in a containers.Map (double->double handles).
ks = m.keys;
for i = 1:numel(ks)
    try Screen('Close', m(ks{i})); end %#ok<TRYNC>
end
end

function cleanUpVideo(win)
    sca;
end

function plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, ~, vbls, targetVBLs, carrierHzs, ifi, areas)
% Aligned with static-image diagnostics: same figures & metrics.

nAreas = size(mod_signals,1);
flickerModeList = string({areas.flickerMode});

% --- Timing vectors (handle NaNs for targetVBLs gracefully) ---
actual_intervals    = diff(vbls);                         % s
scheduled_intervals = diff(targetVBLs);                   % s (may contain NaNs)
timing_error_full   = vbls - targetVBLs;                  % s (may contain NaNs)
valid_te            = isfinite(timing_error_full);
timing_error_valid  = timing_error_full(valid_te);        % s

% --- Summary stats ---
nomHz     = 1/ifi;
effHz     = 1/mean(actual_intervals);
jitter_ms = std(actual_intervals)*1000;
drops     = sum(actual_intervals > 1.5*ifi);
mean_sched_int = mean(scheduled_intervals(isfinite(scheduled_intervals)));
schedHz   = 1/mean_sched_int;

fprintf('\n=== Frame Timing Diagnostics ===\n');
fprintf('Flicker modes & carriers:\n');
for k = 1:nAreas
    chz = carrierHzs(min(k, numel(carrierHzs)));
    fprintf('  Area %d: mode=%s, carrier=%.3f Hz\n', k, flickerModeList(k), chz);
end
fprintf('Mean actual interval: %.5f s (%.2f Hz), SD: %.5f ms\n', ...
    mean(actual_intervals), effHz, std(actual_intervals)*1000);
if isfinite(mean_sched_int)
    fprintf('Mean scheduled interval: %.5f s (%.2f Hz)\n', mean_sched_int, schedHz);
else
    fprintf('Mean scheduled interval: n/a (no finite targets)\n');
end
if ~isempty(timing_error_valid)
    fprintf('Mean abs. timing error: %.5f ms (SD: %.5f ms)\n', ...
        mean(abs(timing_error_valid))*1000, std(timing_error_valid)*1000);
else
    fprintf('Mean abs. timing error: n/a (no finite targets)\n');
end
fprintf('Drops (>1.5×IFI): %d\n', drops);

% --- Figure layout (aligned to static diagnostics) ---
figure('Name','Flicker Diagnostics','NumberTitle','off');
nCols = max(2, nAreas);
tl = tiledlayout(5, nCols, 'TileSpacing','compact');

% Row 1: code sequences
for k = 1:nAreas
    nexttile;
    if ~isempty(code_long_all{k})
        stairs(1:numel(code_long_all{k}), code_long_all{k}, 'LineWidth', 1.1);
        ylim([-0.2 1.2]); ylabel('Code'); xlabel('Frame');
        title(sprintf('Area %d: code', k));
    else
        plot(nan); ylim([0 1]); ylabel('Code'); xlabel('Frame');
        title(sprintf('Area %d: code (none)', k));
    end
    grid on;
end

% Row 2: luminance time series
for k = 1:nAreas
    nexttile;
    plot(t, all_mod_lum(k,:), 'LineWidth', 1.1); grid on;
    fm = flickerModeList(k);
    chz = carrierHzs(min(k, numel(carrierHzs)));
    title(sprintf('Area %d: lum (%s, %.2f Hz)', k, fm, chz));
    xlabel('Time (s)'); ylabel('Lum.');
end

% Row 3: autocorrelation per area
for k = 1:nAreas
    nexttile;
    [acf, lags] = xcorr(mod_signals(k,:), 'coeff');
    plot(lags, acf, 'LineWidth', 1.1); grid on; ylim([0 1]);
    title(sprintf('Area %d: autocorr', k));
    xlabel('Lag (frames)'); ylabel('Norm. corr');
end

% Row 4: Cross-corr (col 1) + Timing-error histogram (col 2..end)
nexttile;
if nAreas >= 2
    [ccf, lagsC] = xcorr(mod_signals(1,:), mod_signals(2,:), 'coeff');
    plot(lagsC, ccf, 'LineWidth', 1.1); grid on;
    title('Cross-corr: Area 1 vs Area 2');
    xlabel('Lag (frames)'); ylabel('Norm. corr');
else
    plot(nan); grid on; title('Cross-corr (needs ≥2 areas)');
end

nexttile([1, nCols-1]);
if ~isempty(timing_error_valid)
    histogram(timing_error_valid*1000, 30);
else
    histogram(0);  % placeholder
end
grid on; xlabel('Timing Error (ms)'); ylabel('Count');
title('Histogram of Timing Error');

% Row 5: Frame intervals (col 1..end-1) + Flip error over time (last col)
nexttile([1, max(1, nCols-1)]);
h1 = plot(actual_intervals*1000,'-o'); hold on;
h2 = plot(scheduled_intervals*1000,'--');
h3 = yline(ifi*1000, 'k-');
grid on; ylabel('Frame Interval (ms)');
legend([h1 h2 h3], {'Actual','Scheduled',sprintf('IFI=%.2f ms',ifi*1000)}, 'Location','best');
title('Frame intervals');

nexttile;
plot(timing_error_full*1000,'-o'); grid on;
ylabel('Timing Error (ms)'); xlabel('Frame #');
title('Flip - Scheduled Time');

% Super title with unified stats
title(tl, sprintf('Diagnostics | nAreas=%d | Nominal=%.2f Hz | Achieved=%.2f Hz | Jitter=%.2f ms | Drops=%d', ...
                  nAreas, nomHz, effHz, jitter_ms, drops));
end


function code_smooth = raised_cosine_smooth(code_long, ramp_len)
% Raised-cosine smoothing across transitions in a binary code sequence.
    code_smooth = code_long;
    N = numel(code_long);
    if ramp_len <= 0, return; end
    w = 0.5 * (1 - cos(linspace(0, pi, ramp_len))); % rising edge
    for i = 2:N
        if code_long(i) ~= code_long(i-1)
            if i-ramp_len+1 < 1, continue; end
            if code_long(i) == 1
                code_smooth(i-ramp_len+1:i) = w;
            else
                code_smooth(i-ramp_len+1:i) = 1 - w;
            end
        end
    end
end
