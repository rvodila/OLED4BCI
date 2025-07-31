function runStimulusExperiment
    % Number of runs/repeats
    k = 5; % <-- Set this as desired

    % --- Set up experiment/session parameters ---
    prm.exp.name        = 'Flicker Tagging Demo';
    prm.exp.datetime    = datetime;
    prm.screen.width_cm = 53;        % Physical width of screen (cm)
    prm.screen.view_dist_cm = 30;    % Viewing distance (cm)
    prm.screen.bg_color = [255 255 255];

    % --- Psychtoolbox Init ---
    Screen('Preference', 'SkipSyncTests', 1); % Remove for real experiments!
    PsychDefaultSetup(2);
    oldPriority = Priority(2);
    screens = Screen('Screens');
    screenNumber = max(screens);
    [win, winRect] = Screen('OpenWindow', screenNumber, prm.screen.bg_color);

    % Get refresh
    ifi = Screen('GetFlipInterval', win);

    % --- Physical <-> Pixel Conversion ---
    win_w_px = winRect(3);
    win_h_px = winRect(4);
    prm.stim.pixel_size = prm.screen.width_cm / win_w_px;
    prm.stim.pix2deg = @(px) (360/pi * atan(px * prm.stim.pixel_size / (2*prm.screen.view_dist_cm)));
    prm.stim.deg2pix = @(deg) (2 * prm.screen.view_dist_cm * tan(pi/360 * deg) / prm.stim.pixel_size);

    % --- Stimulus parameters ---
    stimNames = {'kakadu.png', 'zebra2.png'}; % List your images here
    nStim = numel(stimNames);
    stimSize_px = 400;   % Side length of square (in pixels)
    flickerHzs = [2, 4];   % Flicker frequencies (Hz) for each stimulus
    displayDuration = 10; % seconds
    overlayAlpha = 128;  % Transparency
    lb_lum = 50;   % Minimum luminance (darkest overlay)
    hb_lum = 220;  % Maximum luminance (brightest overlay)

    % --- Stimulus Preparation: Crop, Resize, Normalize, Texture, Position ---
    [textures, dstRects] = prepareStimuli(win, winRect, stimNames, stimSize_px);

    % --- Precompute Tagging/Modulation Signals ---
    totalFrames = round(displayDuration / ifi);
    t = linspace(0, displayDuration, totalFrames); % Time vector for all frames
    tagSigs = zeros(nStim, totalFrames); % [stim, frame]
    for s = 1:nStim
        tagSigs(s,:) = lb_lum + (hb_lum - lb_lum) * (0.5 + 0.5 * sin(2*pi*flickerHzs(s)*t)); % [0,1] sine wave
    end

    Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    try
        for runIdx = 1:k
            % --- Main Flicker Loop ---
            vbl = Screen('Flip', win);
            for frame = 1:totalFrames
                % Background
                Screen('FillRect', win, prm.screen.bg_color);

                % Draw all stimuli and overlays
                for s = 1:nStim
                    % Draw the image
                    Screen('DrawTexture', win, textures(s), [], dstRects(:,s));

                    % Overlay with current luminance tag
                    overlayLum = round(tagSigs(s, frame));
                    overlayColor = [overlayLum overlayLum overlayLum overlayAlpha];
                    Screen('FillRect', win, overlayColor, dstRects(:,s));
                end

                vbl = Screen('Flip', win);

                % ESC to exit
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyIsDown && keyCode(KbName('ESCAPE'))
                    cleanup(win, textures);
                    Priority(oldPriority);
                    return;
                end
            end

            % --- Inter-run screen (except after last run) ---
            if runIdx < k
                Screen('FillRect', win, prm.screen.bg_color);
                message = sprintf('Run %d of %d complete.\n\nNext run is starting.\nPress SPACE to continue.', runIdx, k);
                DrawFormattedText(win, message, 'center', 'center', [0 0 0]);
                Screen('Flip', win);

                % Wait for SPACE to continue or ESC to quit
                while true
                    [keyIsDown, ~, keyCode] = KbCheck;
                    if keyIsDown && keyCode(KbName('SPACE'))
                        break;
                    elseif keyIsDown && keyCode(KbName('ESCAPE'))
                        cleanup(win, textures);
                        Priority(oldPriority);
                        return;
                    end
                end
                KbReleaseWait; % Wait for key release before next run
            end
        end

        cleanup(win, textures);
        Priority(oldPriority);

    catch ME
        cleanup(win, textures);
        Priority(oldPriority);
        rethrow(ME);
    end
end

% ----------------- Helper Functions -----------------

function [textures, dstRects] = prepareStimuli(win, winRect, imgFiles, targetSquareSize)
    % Loads images, normalizes, makes textures, evenly spaces squares horizontally
    nStim = numel(imgFiles);
    textures = zeros(1, nStim);
    dstRects = zeros(4, nStim);
    yPos = winRect(4) / 2;
    xPos = linspace(winRect(3)/(nStim+1), winRect(3)*nStim/(nStim+1), nStim);

    for s = 1:nStim
        img = imread(imgFiles{s});

        % Remove alpha channel if present
        if size(img,3) == 4
            img = img(:,:,1:3);
        end

        % Force grayscale to RGB
        if size(img,3) == 1
            img = repmat(img, [1 1 3]);
        end

        % Center-crop to square
        sz = size(img);
        minDim = min(sz(1), sz(2));
        rowStart = floor((sz(1)-minDim)/2) + 1;
        colStart = floor((sz(2)-minDim)/2) + 1;
        imgSq = img(rowStart:rowStart+minDim-1, colStart:colStart+minDim-1, :);

        % Resize to target square size
        imgSq = imresize(imgSq, [targetSquareSize targetSquareSize]);

        % Force uint8, [0,255]
        imgSq = im2uint8(mat2gray(imgSq));

        % Create PTB texture
        textures(s) = Screen('MakeTexture', win, imgSq);

        % Calculate evenly spaced X positions, centered vertically
        dstRects(:,s) = CenterRectOnPointd([0 0 targetSquareSize targetSquareSize], xPos(s), yPos);
    end
end

function cleanup(win, textures)
    for i = 1:numel(textures)
        Screen('Close', textures(i));
    end
    sca;
end
