% @author: Radovan Vodila (radovan.vodila@ru.nl)


function video_tagging

    %% ---- PARAMETERS ----
    flickerMode = 'hybrid'; % 'freq', 'code', or 'hybrid'
    overlayAlpha = 128;
    lb_lum = 60; hb_lum = 200;
    framesPerBit = 1;
    carrierHzs = [3, 1];
    maxDisplaySec = 5;
    ramp_len = 2;  % param for raised cosine smoothing
    % Stimulus
    rectW = 300; rectH = 150;
    rel_xs = [0.4, 0.6]; rel_ys = [0.3, 0.8];

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
    %HideCursor;

    ifi = Screen('GetFlipInterval', win)
    %ifi = 0.01672; % for testing with 60Hz display setup returnign a non-valid ifi value
    totalFrames = max(round(maxDisplaySec / ifi), 1);
    t = linspace(0, maxDisplaySec, totalFrames);

    %% ---- PRECOMPUTE MODULATION ----
    nOverlays = 2;
    all_mod_lum = zeros(nOverlays, totalFrames);
    codes = {code, code2};
    mod_signals = zeros(nOverlays, totalFrames);
    code_long_all = cell(1, nOverlays);

    for k = 1:nOverlays
        cur_code = codes{k};                              % binary [0,1]
        code_expanded = repelem(cur_code, framesPerBit);
        nrep = ceil(totalFrames / numel(code_expanded));
        code_long = repmat(code_expanded, 1, nrep);
        code_long = code_long(1:totalFrames);

        % --- Raised cosine smoothing (for binary, i.e. [0,1])
        pad_val = code_long(1);
        code_long_padded = [repmat(pad_val, 1, ramp_len) code_long];
        code_long_smoothed = raised_cosine_smooth(code_long_padded, ramp_len);
        code_long = code_long_smoothed(ramp_len+1:end); % remove padding
        code_long_all{k} = code_long;                    % store *binary* smoothed code

        % --- Conversion to bipolar ONLY for modulation
        code_bipolar = 2*code_long - 1;                  % bipolar: [-1,1]

        carrier = 0.5 + 0.5 * sin(2*pi*carrierHzs(k)*t);

        switch lower(flickerMode)
            case 'freq'
                mod_signal = carrier;
            case 'code'
                mod_signal = code_long;                  % still [0,1], pure binary
            case 'hybrid'
                mod_signal = carrier .* code_bipolar;    % now [-1,1], full contrast
            otherwise
                error('Unknown flickerMode: %s', flickerMode);
        end

        mod_signals(k,:) = mod_signal; % store raw mod signal (for diagnostics)
        all_mod_lum(k,:) = lb_lum + (hb_lum - lb_lum) * (mod_signal+1)/2; % [0,1] mapping
    end

    %% ---- VIDEO STIMULUS LOOP ----

    frameCount = 0;
    newFrameCount = 0; % keeps track of actually changed frame content
    dstRect = [];
    vidW = []; vidH = [];

    vbls = zeros(1, totalFrames);
    targetVBLs = zeros(1, totalFrames);

    prevTex = [];

    movie = Screen('OpenMovie', win, movieFile);
    Screen('PlayMovie', movie, 1, [], 1.0);

    % Start timing
    vbl = Screen('Flip', win);
    try
        while true
            tex = Screen('GetMovieImage', win, movie, 0); 
            if tex > 0
                % new frame → close previous texture
                if ~isempty(prevTex)
                    Screen('Close', prevTex);
                end
                prevTex = tex;
                newFrameCount = newFrameCount + 1;
            elseif tex < 0
                % end of movie reached
                break;
            end
                
            if ~isempty(prevTex)            

                frameCount = frameCount + 1;

                % On first frame, get video size
                if isempty(dstRect)
                    rect = Screen('Rect', tex);
                    vidW = rect(3) - rect(1);
                    vidH = rect(4) - rect(2);
                    winW = winRect(3); winH = winRect(4);
                    vidRectW = winW / 2;
                    vidRectH = winH / 2;
                    dstRect = CenterRectOnPointd([0 0 vidRectW vidRectH], winW/2, winH/2);

                    % Overlay rectangles: position relative to displayed video
                    overlayRects = zeros(4,2);
                    for k = 1:2
                        x = dstRect(1) + rel_xs(k)* (dstRect(3)-dstRect(1));
                        y = dstRect(2) + rel_ys(k)* (dstRect(4)-dstRect(2));
                        overlayRects(:,k) = CenterRectOnPointd([0 0 rectW rectH], x, y);
                    end                    
                end
    
                Screen('DrawTexture', win, prevTex, [], dstRect);

                idx = min(newFrameCount, totalFrames);
                for k = 1:2
                    overlayLum = round(all_mod_lum(k, idx));
                    overlayColor = [overlayLum overlayLum overlayLum overlayAlpha];
                    Screen('FillRect', win, overlayColor, overlayRects(:,k));
                end
                
                targetTime = vbl + ifi;
                vbl = Screen('Flip', win, targetTime);

                vbls(frameCount) = vbl;
                targetVBLs(frameCount) = targetTime;
            end

            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE'))
                break;
            end
            if frameCount >= totalFrames
                break;
            end
        end
    catch ME
        cleanUpVideo(win, movie);
        rethrow(ME);
    end

    cleanUpVideo(win, movie);

    %% ---- DIAGNOSTIC PLOTS ----
    plot_modulation_diagnostics(mod_signals, all_mod_lum, t, code_long_all, flickerMode, vbls, targetVBLs, carrierHzs, ifi);
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

    % --- 1. Code sequences used ---
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

    % --- 2. Luminance modulations ---
    subplot(nRows, nCols, sp);
    plot(t, all_mod_lum(1,:), 'r', 'LineWidth', 1.2);
    title(sprintf('Left: Luminance (%.2f Hz)', carrierHzs(1)));
    xlabel('Time (s)'); ylabel('Lum.'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    plot(t, all_mod_lum(2,:), 'b', 'LineWidth', 1.2);
    title(sprintf('Right: Luminance (%.2f Hz)', carrierHzs(2)));
    xlabel('Time (s)'); ylabel('Lum.'); grid on; sp = sp+1;

    % --- 3. Autocorrelations ---
    subplot(nRows, nCols, sp);
    [acfL, lagsL] = xcorr(mod_signals(1,:), 'coeff');
    plot(lagsL, acfL, 'r', 'LineWidth', 1.2);
    title('Left: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on; sp = sp+1;

    subplot(nRows, nCols, sp);
    [acfR, lagsR] = xcorr(mod_signals(2,:), 'coeff');
    plot(lagsR, acfR, 'b', 'LineWidth', 1.2);
    title('Right: Autocorrelation');
    xlabel('Lag (frames)'); ylabel('Norm. corr'); grid on; sp = sp+1;

    % --- 4. Cross-correlation & histogram of timing error ---
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

    % --- 5. Frame interval plot with IFI line ---
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
    % code_long: vector of 0/1 (your upsampled code)
    % ramp_len: number of frames to smooth at each transition
    code_smooth = code_long;
    N = numel(code_long);
    w = 0.5 * (1 - cos(linspace(0, pi, ramp_len))); % rising edge
    for i = 2:N
        if code_long(i) ~= code_long(i-1)
            if code_long(i) == 1  % rising edge: 0 to 1
                code_smooth(i-ramp_len+1:i) = w;
            else                 % falling edge: 1 to 0
                code_smooth(i-ramp_len+1:i) = 1-w;
            end
        end
    end
end
