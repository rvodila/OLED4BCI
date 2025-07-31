"""
@author: Radovan Vodila (radovan.vodila@ru.nl)
"""

function flicker_protocol_video_hybrid

    %% ---- USER OPTIONS ----
    flickerMode = 'hybrid'; % 'freq', 'code', or 'hybrid'
    overlayAlpha = 128;
    lb_lum = 0; hb_lum = 255;
    bitPerFrame = 1;
    carrierHzs = [24, 25];
    maxDisplaySec = 5;

    % --- Video file ---
    movieFile = fullfile(pwd, 'project', 'stimulus', 'images', 'ape_walk.mp4');

    % --- Code files ---
    codefile = fullfile(pwd, 'project', 'stimulus', 'files', 'mgold_61_6521.mat');
    S = load(codefile);
    code  = double(S.codes(1, :));
    code2 = double(S.codes(2, :));
    code  = code(:)'; code2 = code2(:)'; % row vectors

    %% ---- PSYCHTOOLBOX SETUP ----
    Screen('Preference', 'SkipSyncTests', 1);
    PsychDefaultSetup(2);
    screens = Screen('Screens');
    screenNumber = max(screens);
    bgColor = [128 128 128];
    [win, winRect] = Screen('OpenWindow', screenNumber, bgColor);
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;

    ifi = Screen('GetFlipInterval', win);
    totalFrames = max(round(maxDisplaySec / ifi), 1);
    t = linspace(0, maxDisplaySec, totalFrames);

    % ---- PRECOMPUTE MODULATION ----
    nOverlays = 2;
    all_mod_lum = zeros(nOverlays, totalFrames);
    codes = {code, code2};
    mod_signals = zeros(nOverlays, totalFrames);

    for k = 1:nOverlays
        cur_code = codes{k};
        code_expanded = repelem(cur_code, bitPerFrame);
        nrep = ceil(totalFrames / numel(code_expanded));
        code_long = repmat(code_expanded, 1, nrep);
        code_long = code_long(1:totalFrames);

        % Bipolar mapping: 0 -> -1, 1 -> 1 (for true Gold code properties)
        code_long = 2 * code_long - 1;

        carrier = 0.5 + 0.5 * sin(2*pi*carrierHzs(k)*t);

        switch lower(flickerMode)
            case 'freq'
                mod_signal = carrier;
            case 'code'
                mod_signal = code_long;
            case 'hybrid'
                mod_signal = carrier .* code_long;
            otherwise
                error('Unknown flickerMode: %s', flickerMode);
        end

        mod_signals(k,:) = mod_signal; % store raw mod signal (for diagnostics)
        all_mod_lum(k,:) = lb_lum + (hb_lum - lb_lum) * (mod_signal+1)/2; % map [-1,1] -> [0,1]
    end

    %% ---- VIDEO STIMULUS LOOP ----
    movie = Screen('OpenMovie', win, movieFile);
    rectW = 500; rectH = 300;
    rel_xs = [0.4, 0.6]; rel_ys = [0.3, 0.8];
    Screen('PlayMovie', movie, 1);
    vbl = Screen('Flip', win);
    frameCount = 0;
    dstRect = [];
    vidW = []; vidH = [];

    try
        while true
            tex = Screen('GetMovieImage', win, movie, 1);
            if tex <= 0
                break;
            end
            frameCount = frameCount + 1;

            % On first frame, get video size (fast method)
            if isempty(dstRect)
                rect = Screen('Rect', tex);
                vidW = rect(3) - rect(1);
                vidH = rect(4) - rect(2);
                winW = winRect(3); winH = winRect(4);
                vidRectW = winW / 2;
                vidRectH = winH / 2;
                dstRect = CenterRectOnPointd([0 0 vidRectW vidRectH], winW/2, winH/2);
            end

            Screen('DrawTexture', win, tex, [], dstRect);

            % Overlay rectangles: position relative to displayed video
            overlayRects = zeros(4,2);
            for k = 1:2
                x = dstRect(1) + rel_xs(k)* (dstRect(3)-dstRect(1));
                y = dstRect(2) + rel_ys(k)* (dstRect(4)-dstRect(2));
                overlayRects(:,k) = CenterRectOnPointd([0 0 rectW rectH], x, y);
            end

            idx = min(frameCount, totalFrames);
            for k = 1:2
                overlayLum = round(all_mod_lum(k, idx));
                overlayColor = [overlayLum overlayLum overlayLum overlayAlpha];
                Screen('FillRect', win, overlayColor, overlayRects(:,k));
            end

            vbl = Screen('Flip', win, vbl + 0.5*ifi);
            Screen('Close', tex);

            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE'))
                break;
            end
        end
    catch ME
        cleanUpVideo(win, movie);
        rethrow(ME);
    end

    cleanUpVideo(win, movie);

    %% ---- DIAGNOSTIC PLOTS ----
    plot_modulation_diagnostics(mod_signals, all_mod_lum, t, codes, flickerMode);
end

function cleanUpVideo(win, movie)
    sca;
    try
        if ~isempty(movie) && isnumeric(movie) && movie > 0
            Screen('CloseMovie', movie);
        end
    end
    ShowCursor;
end

function plot_modulation_diagnostics(mod_signals, all_mod_lum, t, codes, flickerMode)
    % Diagnostic plotting for code, luminance, auto/cross-correlation

    figure('Name','Image Flicker Diagnostics', 'NumberTitle', 'off');
    nRows = 4; nCols = 2;
    sp = 1;

    % --- 1. Plot code sequences ---
    subplot(nRows, nCols, sp);
    stairs(1:numel(codes{1}), codes{1}, 'r', 'LineWidth', 1.2);
    ylim([-1.2 1.2]);
    title('Left code sequence');
    ylabel('Bit'); xlabel('Code bit');
    grid on;
    sp = sp + 1;

    subplot(nRows, nCols, sp);
    stairs(1:numel(codes{2}), codes{2}, 'b', 'LineWidth', 1.2);
    ylim([-1.2 1.2]);
    title('Right code sequence');
    ylabel('Bit'); xlabel('Code bit');
    grid on;
    sp = sp + 1;

    % --- 2. Plot luminance modulations ---
    subplot(nRows, nCols, sp);
    plot(t, all_mod_lum(1,:), 'r', 'LineWidth', 1.2);
    title('Left: Luminance Modulation');
    xlabel('Time (s)'); ylabel('Lum.');
    grid on;
    sp = sp + 1;

    subplot(nRows, nCols, sp);
    plot(t, all_mod_lum(2,:), 'b', 'LineWidth', 1.2);
    title('Right: Luminance Modulation');
    xlabel('Time (s)'); ylabel('Lum.');
    grid on;
    sp = sp + 1;

    % --- 3. Plot autocorrelations ---
    subplot(nRows, nCols, sp);
    [acfL, lagsL] = xcorr(mod_signals(1,:), 'coeff');
    plot(lagsL, acfL, 'r', 'LineWidth', 1.2);
    title('Left: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr');
    grid on;
    sp = sp + 1;

    subplot(nRows, nCols, sp);
    [acfR, lagsR] = xcorr(mod_signals(2,:), 'coeff');
    plot(lagsR, acfR, 'b', 'LineWidth', 1.2);
    title('Right: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr');
    grid on;
    sp = sp + 1;

    % --- 4. Plot cross-correlation (spans both columns) ---
    subplot(nRows, nCols, sp:(nRows*nCols));
    [ccf, lagsC] = xcorr(mod_signals(1,:), mod_signals(2,:), 'coeff');
    plot(lagsC, ccf, 'k', 'LineWidth', 1.2);
    title('Cross-correlation: Left vs Right');
    xlabel('Lag (frames)'); ylabel('Norm. corr');
    grid on;

    % --- Overall title ---
    sgtitle(['Video Flicker Diagnostics: ', flickerMode]);
end

