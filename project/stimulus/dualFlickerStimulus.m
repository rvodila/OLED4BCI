"""
@author: Radovan Vodila (radovan.vodila@ru.nl)
"""

function flicker_protocol_two_images_hybrid
    %% ---- USER OPTIONS ----
    flickerMode = 'hybrid'; % 'freq', 'code', or 'hybrid'
    maxDisplaySec = 5;
    bitPerFrame = 2;
    overlayAlpha = 128;
    lb_lum = 0; hb_lum = 255;
    stimSize = 400;
    rel_xs = [0.4, 0.6]; rel_y = 0.5;
    carrierHzs = [2, 4];
    imageFiles = {fullfile(pwd, 'project', 'stimulus', 'images', 'kakadu.png'), ...
                  fullfile(pwd, 'project', 'stimulus', 'images', 'zebra2.png')};

    %% ---- CODE FILE ----
    codefile = fullfile(pwd, 'project', 'stimulus', 'files', 'mgold_61_6521.mat');
    S = load(codefile);
    code  = double(S.codes(1, :));
    code2 = double(S.codes(2, :));
    code  = code(:)'; code2 = code2(:)';

    %% ---- PSYCHTOOLBOX SETUP ----
    Screen('Preference', 'SkipSyncTests', 1);
    PsychDefaultSetup(2);
    screens = Screen('Screens');
    screenNumber = max(screens);
    bgColor = [255 255 255];
    [win, winRect] = Screen('OpenWindow', screenNumber, bgColor);
    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;

    ifi = Screen('GetFlipInterval', win);
    totalFrames = max(round(maxDisplaySec / ifi), 1);
    t = linspace(0, maxDisplaySec, totalFrames);

    %% ---- LOAD & CENTER IMAGES ----
    nStim = 2;
    textures = zeros(1, nStim);
    dstRects = zeros(4, nStim);
    winW = winRect(3); winH = winRect(4);

    for k = 1:nStim
        img = imread(imageFiles{k});
        if size(img,3)==1
            img = repmat(img, [1 1 3]);
        end
        sz = size(img);
        minDim = min(sz(1:2));
        rowStart = floor((sz(1)-minDim)/2)+1;
        colStart = floor((sz(2)-minDim)/2)+1;
        imgSq = img(rowStart:rowStart+minDim-1, colStart:colStart+minDim-1, :);
        imgSq = imresize(imgSq, [stimSize stimSize]);
        imgSq = im2uint8(mat2gray(imgSq));
        textures(k) = Screen('MakeTexture', win, imgSq);
        x = winW * rel_xs(k); y = winH * rel_y;
        dstRects(:,k) = CenterRectOnPointd([0 0 stimSize stimSize], x, y);
    end

    %% ---- PRECOMPUTE FLICKER MODULATION ----
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
        code_long = 2 * code_long - 1; % Bipolar mapping

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

        mod_signals(k,:) = mod_signal;
        all_mod_lum(k,:) = lb_lum + (hb_lum - lb_lum) * (mod_signal+1)/2; % map [-1,1] -> [0,1]; autocor attributes
    end

    %% ---- STIMULUS LOOP ----
    vbl = Screen('Flip', win);
    frameCount = 0;
    try
        while true
            frameCount = frameCount + 1;
            idx = min(frameCount, totalFrames);

            Screen('FillRect', win, bgColor);
            for k = 1:nStim
                Screen('DrawTexture', win, textures(k), [], dstRects(:,k));
                overlayLum = round(all_mod_lum(k, idx));
                overlayColor = [overlayLum overlayLum overlayLum overlayAlpha];
                Screen('FillRect', win, overlayColor, dstRects(:,k));
            end

            vbl = Screen('Flip', win, vbl + 0.5*ifi);

            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE'))
                break;
            end
            if frameCount >= totalFrames
                break;
            end
        end
    catch ME
        cleanup(win, textures);
        rethrow(ME);
    end
    cleanup(win, textures);

    %% ---- DIAGNOSTIC PLOTS ----
    plot_modulation_diagnostics(mod_signals, all_mod_lum, t, codes, flickerMode);
end

function cleanup(win, textures)
    for i = 1:numel(textures)
        Screen('Close', textures(i));
    end
    sca;
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
    sgtitle(['Image Flicker Diagnostics: ', flickerMode]);
end
