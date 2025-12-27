%% === mainSinFlickerTest.m ===

% --- Setup Psychtoolbox ---
Screen('Preference','SkipSyncTests', 1); % skip for dev
PsychDefaultSetup(2);
KbName('UnifyKeyNames');

screenNumber = max(Screen('Screens'));
bg = 0.5 * 255;  % mid-gray
[window, winRect] = Screen('OpenWindow', screenNumber, bg);

ifi = Screen('GetFlipInterval', window); % measured frame interval

% --- Stimulus parameters ---
numCycles = 50;          % number of 4-frame cycles
DURATION_TRIGGER = NaN;  % no triggers (skip)
address = [];            % dummy argument

% --- Run your flicker function ---
sinusoidalFlickerFinal_240Hz_olaf(window, ifi, numCycles, DURATION_TRIGGER, address);

% --- Clean up ---
sca;
ShowCursor;
Priority(0);
