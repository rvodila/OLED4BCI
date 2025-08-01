"""
@author: Radovan Vodila (radovan.vodila@ru.nl)
"""

function flicker_protocol_two_images_hybrid
    %% ---- USER OPTIONS ----
    flickerMode = 'hybrid'; % 'freq', 'code', or 'hybrid'
    maxDisplaySec = 10;
    framesPerBit = 4;
    overlayAlpha = 128;
    lb_lum = 50; hb_lum = 195;
    stimSize = 400;
    rel_xs = [0.4, 0.6]; rel_y = 0.5;
    carrierHzs = [5, 15];
    imageFiles = {fullfile(pwd, 'project', 'stimulus', 'images', 'kakadu.png'), ...
                  fullfile(pwd, 'project', 'stimulus', 'images', 'zebra2.png')};
    ramp_len = 4; % number of frames over which to smooth transitions

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
    code_long_all = cell(1, nOverlays);
    for k = 1:nOverlays
        cur_code = codes{k};
        code_expanded = repelem(cur_code, framesPerBit);
        nrep = ceil(totalFrames / numel(code_expanded));
        code_long = repmat(code_expanded, 1, nrep);
        code_long = code_long(1:totalFrames);
        % code_long = 2 * code_long - 1; % bipolar

        % --- Padding for smoothing at start
        pad_val = code_long(1);
        code_long_padded = [repmat(pad_val, 1, ramp_len) code_long];
        code_long_smoothed = raised_cosine_smooth(code_long_padded, ramp_len);
        code_long = code_long_smoothed(ramp_len+1:end); % remove padding
        code_long_all{k} = code_long;

        % CARRIER
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
    vbls = zeros(1, totalFrames);
    targetVBLs = zeros(1, totalFrames);

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

                targetTime = vbl + 0.5*ifi;      % What you asked for (scheduled time)
                vbl = Screen('Flip', win, targetTime);   % When the frame actually flipped
                vbls(frameCount) = vbl;                 % Store actual flip time
                targetVBLs(frameCount) = targetTime;    % Store target flip time


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
    % After the stimulus loop, call:
    plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, flickerMode, vbls, targetVBLs, carrierHzs, ifi);



end

function cleanup(win, textures)
    for i = 1:numel(textures)
        Screen('Close', textures(i));
    end
    sca;
    ShowCursor;
end

function plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, flickerMode, vbls, targetVBLs, carrierHzs, ifi)
    % --- Frame timing analysis ---
    actual_intervals = diff(vbls);         % Actual frame durations (s)
    scheduled_intervals = diff(targetVBLs);% What you requested (s)
    timing_error = vbls - targetVBLs;      % Flip error (s)

    fprintf('\n=== Frame Timing Diagnostics ===\n');
    fprintf('Flicker frequencies (Hz): Left: %.3f, Right: %.3f\n', carrierHzs(1), carrierHzs(2));
    fprintf('Mean actual interval: %.5f s (%.2f Hz), SD: %.5f ms\n', ...
        mean(actual_intervals), 1/mean(actual_intervals), std(actual_intervals)*1000);
    fprintf('Mean scheduled interval: %.5f s (%.2f Hz)\n', ...
        mean(scheduled_intervals), 1/mean(scheduled_intervals));
    fprintf('Mean abs. timing error: %.5f ms (SD: %.5f ms)\n', ...
        mean(abs(timing_error))*1000, std(timing_error)*1000);

    nRows = 5; nCols = 2; sp = 1;

    figure('Name','Image Flicker Diagnostics','NumberTitle','off');

    % --- 1. Code sequences ---
% --- 1. Code sequences actually used for modulation ---
    subplot(nRows, nCols, sp);
    stairs(1:numel(code_long_all{1}), code_long_all{1}, 'r', 'LineWidth', 1.2);
    ylim([-1.2 1.2]);
    title('Left: modified code sequence');
    ylabel('Code value'); xlabel('Frame'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    stairs(1:numel(code_long_all{2}), code_long_all{2}, 'b', 'LineWidth', 1.2);
    ylim([-1.2 1.2]);
    title('Right: modified code sequence');
    ylabel('Code value'); xlabel('Frame'); grid on; sp = sp+1;

    % --- 2. Luminance time series ---
    subplot(nRows,nCols,sp);
    plot(t, all_mod_lum(1,:), 'r', 'LineWidth', 1.2);
    title(sprintf('Left: Luminance (%.2f Hz)', carrierHzs(1)));
    xlabel('Time (s)'); ylabel('Lum.'); grid on; sp = sp+1;

    subplot(nRows,nCols,sp);
    plot(t, all_mod_lum(2,:), 'b', 'LineWidth', 1.2);
    title(sprintf('Right: Luminance (%.2f Hz)', carrierHzs(2)));
    xlabel('Time (s)'); ylabel('Lum.'); grid on; sp = sp+1;

    % --- 3. Autocorrelations ---
    subplot(nRows,nCols,sp);
    [acfL, lagsL] = xcorr(mod_signals(1,:), 'coeff');
    plot(lagsL, acfL, 'r', 'LineWidth', 1.2);
    title('Left: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on; sp = sp+1;

    subplot(nRows,nCols,sp);
    [acfR, lagsR] = xcorr(mod_signals(2,:), 'coeff');
    plot(lagsR, acfR, 'b', 'LineWidth', 1.2);
    title('Right: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on; sp = sp+1;

    % --- 4. Cross-correlation (left col), Histogram of timing error (right col) ---
    subplot(nRows,nCols,sp);
    [ccf, lagsC] = xcorr(mod_signals(1,:), mod_signals(2,:), 'coeff');
    plot(lagsC, ccf, 'k', 'LineWidth', 1.2);
    title('Cross-corr: Left vs Right');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on;

    subplot(nRows,nCols,sp+1);
    histogram(timing_error*1000, 30, 'FaceColor', [0.2 0.2 0.8]);
    xlabel('Timing Error (ms)');
    ylabel('Count');
    title('Histogram of Timing Error');

    % --- 5. 
       % --- 5. Frame interval plot ---
    subplot(nRows,nCols,nRows*nCols-1);
    h1 = plot(actual_intervals*1000,'-o'); hold on;
    h2 = plot(scheduled_intervals*1000,'--');
    % Add horizontal line for IFI
    h3 = yline(ifi*1000, 'k-', ...
            'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
    ylabel('Frame Interval (ms)');
    legend([h1 h2 h3], {'Actual', 'Scheduled', sprintf('IFI = %.2f ms', ifi*1000)});
    grid on;
    title('Frame Intervals');

    subplot(nRows,nCols,nRows*nCols);
    plot(timing_error*1000,'-o');
    ylabel('Timing Error (ms)'); xlabel('Frame #');
    title('Flip - Scheduled Time'); grid on;

    sgtitle(['Image Flicker Diagnostics: ', flickerMode]);
end

function code_smooth = raised_cosine_smooth(code_long, ramp_len)
    % code_long: vector of -1/1 (your upsampled code)
    % ramp_len: number of frames to smooth at each transition
    code_smooth = code_long;
    N = numel(code_long);
    w = 0.5 * (1 - cos(linspace(0, pi, ramp_len))); % rising edge
    for i = 2:N
        if code_long(i) ~= code_long(i-1)
            if code_long(i) == 1  % rising edge: -1 to 1
                code_smooth(i-ramp_len+1:i) = w;
            else                 % falling edge: 1 to -1
                code_smooth(i-ramp_len+1:i) = 1-w;
                % code_smooth(i-ramp_len+1:i) = 2*code_smooth(i-ramp_len+1:i)-1; % map to [-1,1], BIPOLAR
            end
        end
    end
end

