function simple_multi_flicker(squares, durationSec)
% SIMPLE_MULTI_FLICKER
% --------------------
% Presents one or more squares, each flickering sinusoidally at its own
% temporal frequency. All drawn per refresh, VBL-locked.
%
% Usage:
%   simple_multi_flicker();                 % default 2 squares
%
%   % custom: two squares, 10s total
%   sq(1) = struct('freq', 12, 'x', 600, 'y', 540, ...
%                  'size', 300, 'lb', 40, 'hb', 200);
%   sq(2) = struct('freq', 15, 'x', 1100, 'y', 540, ...
%                  'size', 300, 'lb', 60, 'hb', 220);
%   simple_multi_flicker(sq, 10);
%
% Press ESC to quit early.

    % -------------------- defaults --------------------
    if nargin < 1 || isempty(squares)
        % Two default squares as def if none given
        squares(1) = struct('freq', 5, 'x', 640-200, 'y', 540, ...
                            'size', 50, 'lb', 40, 'hb', 200);
        squares(2) = struct('freq', 10, 'x', 640-100, 'y', 540, ...
                            'size', 50, 'lb', 40, 'hb', 200);
        squares(3) = struct('freq', 15, 'x', 640+200, 'y', 540, ...
                            'size', 50, 'lb', 40, 'hb', 200);
        squares(4) = struct('freq', 20, 'x', 640+100, 'y', 540, ...
                            'size', 50, 'lb', 40, 'hb', 200);
    end
    if nargin < 2 || isempty(durationSec)
        durationSec = 10;
    end

    devModeSkipSync = false;  % skip sync for dev
    if devModeSkipSync
        Screen('Preference','SkipSyncTests',1);
    else
        Screen('Preference','SkipSyncTests',0);
    end

    PsychDefaultSetup(2); KbName('UnifyKeyNames');
    bg = 128;
    screenNumber = max(Screen('Screens'));
    [win, winRect] = Screen('OpenWindow', screenNumber, bg);
    HideCursor;
    Screen('BlendFunction', win, 'GL_ONE', 'GL_ZERO'); % no blending

    ifi = Screen('GetFlipInterval', win);
    fps = 1/ifi;
    totalFrames = round(durationSec * fps);
    fprintf('%.2f Hz refresh, %d frames total (~%.2f s)\n', fps, totalFrames, totalFrames*ifi);

    % -------------------- build rects & timing --------------------
    nSq = numel(squares);
    baseRects = zeros(4, nSq);
    dstRects  = zeros(4, nSq);

    [cx, cy] = RectCenter(winRect);
    for k = 1:nSq
        if ~isfield(squares(k),'x'), squares(k).x = cx; end
        if ~isfield(squares(k),'y'), squares(k).y = cy; end
        if ~isfield(squares(k),'size'), squares(k).size = 300; end
        if ~isfield(squares(k),'lb'), squares(k).lb = 40; end
        if ~isfield(squares(k),'hb'), squares(k).hb = 200; end
        if ~isfield(squares(k),'freq'), squares(k).freq = 10; end

        baseRects(:,k) = [0 0 squares(k).size squares(k).size];
        dstRects(:,k)  = CenterRectOnPointd(baseRects(:,k), squares(k).x, squares(k).y);
    end

    % time vector at refresh samples
    t = (0:totalFrames-1) * ifi;

    % precompute luminances for each square
    lum = zeros(nSq, totalFrames);
    for k = 1:nSq
        sq = squares(k);
        lum(k,:) = sq.lb + (sq.hb - sq.lb) * (0.5 + 0.5 * sin(2*pi*sq.freq*t));
    end

    % -------------------- stimulus loop --------------------
    Priority(MaxPriority(win));
    vbl = Screen('Flip', win);
    waitframes = 1;

    try
        for f = 1:totalFrames
            Screen('FillRect', win, bg);

            % draw each square with its luminance this frame
            for k = 1:nSq
                rgb = repmat(uint8(round(lum(k,f))), 3, 1);
                Screen('FillRect', win, rgb, dstRects(:,k));
            end

            vbl = Screen('Flip', win, vbl + (waitframes - 0.5)*ifi);

            % quit with ESC
            [down,~,kc] = KbCheck;
            if down && kc(KbName('ESCAPE')), break; end
        end
    catch ME
        sca; ShowCursor; Priority(0);
        rethrow(ME);
    end

    % -------------------- cleanup --------------------
    sca; ShowCursor; Priority(0);
    fprintf('Done.\n');
end
