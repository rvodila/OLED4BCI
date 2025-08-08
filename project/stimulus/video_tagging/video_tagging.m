% @author: Radovan Vodila (radovan.vodila@ru.nl)
function flicker_protocol_video_hybrid

    %% ---- PARAMETERS ----
    flickerMode = 'hybrid'; % 'freq', 'code', or 'hybrid'
    overlayAlpha = 128;
    lb_lum = 60; hb_lum = 200;
    framesPerBit = 1;
    carrierHzs = [3, 1];
    maxDisplayLen = 10;
    ramp_len = 2;  % for raised cosine smoothing
    rectW = 300; rectH = 150;
    rel_xs = [0.4, 0.6]; rel_ys = [0.3, 0.8];
    movieFile = fullfile(pwd, 'project', 'stimulus', 'images', 'ape_walk.mp4'); % 25hz, 17sec, 950 x 540
    codefile = fullfile(pwd, 'project', 'stimulus', 'codes', 'mgold_61_6521.mat');
    S = load(codefile);
    code  = double(S.codes(1, :)); code2 = double(S.codes(2, :));
    code  = code(:)'; code2 = code2(:)';

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
    totalFrames = max(round(maxDisplayLen / ifi), 1);
    t = linspace(0, maxDisplayLen, totalFrames);
    
    %% ---- PRECOMPUTE MODULATION ----
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

        pad_val = code_long(1);
        code_long_padded = [repmat(pad_val, 1, ramp_len) code_long];
        code_long_smoothed = raised_cosine_smooth(code_long_padded, ramp_len);
        code_long = code_long_smoothed(ramp_len+1:end); % remove padding
        code_long_all{k} = code_long;

        code_bipolar = 2*code_long - 1;  % for hybrid

        carrier = 0.5 + 0.5 * sin(2*pi*carrierHzs(k)*t);

        switch lower(flickerMode)
            case 'freq'
                mod_signal = carrier;
            case 'code'
                mod_signal = code_long;         % [0,1]
            case 'hybrid'
                mod_signal = carrier .* code_bipolar; % [-1,1]
            otherwise
                error('Unknown flickerMode: %s', flickerMode);
        end

        mod_signals(k,:) = mod_signal;
        all_mod_lum(k,:) = lb_lum + (hb_lum - lb_lum) * (mod_signal+1)/2;
    end

    %% ---- PRELOAD VIDEO TO RAM (as Textures) ----

    % --- VIDEO SETUP ---
    videoReader = VideoReader(movieFile);
    videoFPS = videoReader.FrameRate;
    videoDuration = videoReader.Duration;
    nVidFrames = floor(videoFPS * videoDuration);
    % --- PTB DISPLAY SETUP ---
    ifi = Screen('GetFlipInterval', win);
    displayFPS = 1/ifi;

    % How many video frames fit in the allotted display time?
    nDisplayableFrames = min(round(videoFPS * maxDisplayLen), nVidFrames);

    % Corresponding display duration
    actualDisplayLen = nDisplayableFrames / videoFPS;

    % Compute number of screen refreshes we will use (best: sync to video FPS)
    % You can either show each video frame for N display frames, or show at real time

    % If you want to match video timing (recommended):
    videoFrameTimes = (0:nDisplayableFrames-1) / videoFPS; % seconds

    % For each display frame (at your display FPS), find which video frame to show
    displayTimes = (0:round(actualDisplayLen*displayFPS)-1) / displayFPS;

    % For each display frame, which video frame to show?
    videoFrameIdx = interp1(videoFrameTimes, 1:nDisplayableFrames, displayTimes, 'nearest', 'extrap');
    videoFrameIdx = min(max(round(videoFrameIdx),1),nVidFrames); % safety

    % --- LOAD VIDEO TO TEXTURES
    videoFrames = cell(1, nVidFrames);
    for f = 1:nVidFrames
        frame = readFrame(videoReader);
        if size(frame,3)==1 % if single channel, convert to RGB
            frame = repmat(frame, [1 1 3]);
        end
        videoFrames{f} = im2uint8(frame);
    end

    % Print memory/resource usage:
    infoFrames = whos('videoFrames');
    fprintf('\n[Resource Usage] videoFrames: %.2f MB (%.2f GB)\n', ...
        infoFrames.bytes/2^20, infoFrames.bytes/2^30);
    vars = whos; totalBytes = sum([vars.bytes]);
    fprintf('[Resource Usage] Total workspace: %.2f MB (%.2f GB)\n', ...
        totalBytes/2^20, totalBytes/2^30);
    mem = memory;
    fprintf('[System RAM] Used: %.2f GB | Free: %.2f GB | Total: %.2f GB\n', ...
    mem.MemUsedMATLAB/2^30, mem.MemAvailableAllArrays/2^30, mem.MaxPossibleArrayBytes/2^30);

    % Convert videoFrames to Psychtoolbox textures (pre-make all!)
    videoTextures = zeros(1, nVidFrames);
    for f = 1:nVidFrames
        videoTextures(f) = Screen('MakeTexture', win, videoFrames{f});
    end
    clear videoFrames; % save RAM after making textures

    %% ---- VIDEO STIMULUS LOOP (Precomputed Textures) ----
    frameCount = 0;
    vbl = Screen('Flip', win);
    vbls = zeros(1, totalFrames);
    targetVBLs = zeros(1, totalFrames);

    winW = winRect(3); winH = winRect(4);
    vidRectW = winW / 2;
    vidRectH = winH / 2;
    dstRect = CenterRectOnPointd([0 0 vidRectW vidRectH], winW/2, winH/2);
    try
        for frameCount = 1:totalFrames
            vidIdx = videoFrameIdx(frameCount);
            if vidIdx > nVidFrames
                break;
            end 
            tex = videoTextures(vidIdx);
            Screen('DrawTexture', win, tex, [], dstRect);

            overlayRects = zeros(4,2);
            for k = 1:2
                x = dstRect(1) + rel_xs(k)* (dstRect(3)-dstRect(1));
                y = dstRect(2) + rel_ys(k)* (dstRect(4)-dstRect(2));
                overlayRects(:,k) = CenterRectOnPointd([0 0 rectW rectH], x, y);
            end

            for k = 1:2
                overlayLum = round(all_mod_lum(k, frameCount));
                overlayColor = [overlayLum overlayLum overlayLum overlayAlpha];
                Screen('FillRect', win, overlayColor, overlayRects(:,k));
            end

            targetTime = vbl + 0.5*ifi;
            vbl = Screen('Flip', win, targetTime);
            vbls(frameCount) = vbl;
            targetVBLs(frameCount) = targetTime;

            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE'))
                break;
            end
            if frameCount >= totalFrames
                break;
            end
        end
        Screen('Close', videoTextures); % This closes all in one call!
        cleanUpVideo(win);
    catch ME
        Screen('Close', videoTextures); % Safe to call even if already closed, just does nothing
        cleanUpVideo(win);
        rethrow(ME);
    end

    %% ---- DIAGNOSTIC PLOTS ----
    plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, flickerMode, vbls, targetVBLs, carrierHzs, ifi);
end

function cleanUpVideo(win)
    sca;
    ShowCursor;
end

function plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, flickerMode, vbls, targetVBLs, carrierHzs, ifi)
    actual_intervals = diff(vbls);
    scheduled_intervals = diff(targetVBLs);
    timing_error = vbls - targetVBLs;
    fprintf('\n=== Frame Timing Diagnostics ===\n');
    fprintf('Flicker frequencies (Hz): Left: %.3f, Right: %.3f\n', carrierHzs(1), carrierHzs(2));
    fprintf('Mean actual interval: %.5f s (%.2f Hz), SD: %.5f ms\n', ...
        mean(actual_intervals), 1/mean(actual_intervals), std(actual_intervals)*1000);
    fprintf('Mean scheduled interval: %.5f s (%.2f Hz)\n', ...
        mean(scheduled_intervals), 1/mean(scheduled_intervals));
    fprintf('Mean abs. timing error: %.5f ms (SD: %.5f ms)\n', ...
        mean(abs(timing_error))*1000, std(timing_error)*1000);

    nRows = 5; nCols = 2; sp = 1;
    figure('Name','Video Flicker Diagnostics','NumberTitle','off');

    subplot(nRows, nCols, sp);
    stairs(1:numel(code_long_all{1}), code_long_all{1}, 'r', 'LineWidth', 1.2);
    ylim([-0.2 1.2]);
    title('Left: code sequence used');
    ylabel('Code value'); xlabel('Frame'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    stairs(1:numel(code_long_all{2}), code_long_all{2}, 'b', 'LineWidth', 1.2);
    ylim([-0.2 1.2]);
    title('Right: code sequence used');
    ylabel('Code value'); xlabel('Frame'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    plot(t, all_mod_lum(1,:), 'r', 'LineWidth', 1.2);
    title(sprintf('Left: Modulated luminance (%.2f Hz)', carrierHzs(1)));
    xlabel('Time (s)'); ylabel('Lum.'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    plot(t, all_mod_lum(2,:), 'b', 'LineWidth', 1.2);
    title(sprintf('Right: Modulated Luminance (%.2f Hz)', carrierHzs(2)));
    xlabel('Time (s)'); ylabel('Lum.'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    [acfL, lagsL] = xcorr(mod_signals(1,:), 'coeff');
    plot(lagsL, acfL, 'r', 'LineWidth', 1.2);
    ylim([0 1]);
    title('Left: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    [acfR, lagsR] = xcorr(mod_signals(2,:), 'coeff');
    plot(lagsR, acfR, 'b', 'LineWidth', 1.2);
    ylim([0 1]);
    title('Right: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on; sp = sp+1;

    subplot(nRows,nCols,sp);
    [ccf, lagsC] = xcorr(mod_signals(1,:), mod_signals(2,:), 'coeff');
    plot(lagsC, ccf, 'k', 'LineWidth', 1.2);
    ylim([0 1]);
    title('Cross-corr: Left vs Right');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on;

    subplot(nRows,nCols,sp+1);
    histogram(timing_error*1000, 30, 'FaceColor', [0.2 0.2 0.8]);
    xlabel('Timing Error (ms)');
    ylabel('Count');
    title('Histogram of Timing Error');

    subplot(nRows,nCols,nRows*nCols-1);
    h1 = plot(actual_intervals*1000,'-o'); hold on;
    h2 = plot(scheduled_intervals*1000,'--');
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

    sgtitle(['Video Flicker Diagnostics: ', flickerMode]);
end

function code_smooth = raised_cosine_smooth(code_long, ramp_len)
    code_smooth = code_long;
    N = numel(code_long);
    w = 0.5 * (1 - cos(linspace(0, pi, ramp_len))); % rising edge
    for i = 2:N
        if code_long(i) ~= code_long(i-1)
            if code_long(i) == 1
                code_smooth(i-ramp_len+1:i) = w;
            else
                code_smooth(i-ramp_len+1:i) = 1-w;
            end
        end
    end
end
