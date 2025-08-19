Dfunction sinusoidalFlickerFinal_240Hz(window, ifi, numCycles, DURATION_TRIGGER, address)

%% (Pseudo)Sinusoidal flicker for testing 240 Hz OLED response times
% This function creates a pseudosinusoidal flicker that iterates from black
% to white and back to black with an intermittent grey level, the flicker
% is tapered at the beginning and end of the desired presentation sequence
% (length determined by the function input numCycles) according to a 
% Hann/raised cosine taper. Written by Arne STEIN & Olaf DIMIGEN, 2024

%% Create Hanning Taper 
% Here we create a Hanning taper for the sequence
waitframes = 1;

hannTaper = hann(96)'; % Hanning taper for 240hz screen 
seqLum = repmat([-100 0 100 0], 1, 24); % 24 repetitions because 96 taper points
finTaper = seqLum .* hannTaper; % Multiply contrast values with taper

% Determine stimulus size
si = 64;

% Size of support in pixels, derived from si:
tw = 20*si+1;
th = 20*si+1;

backgroundcolorOffset = [0.5 0.5 0.5 0]; % Background grey level
blobtex = CreateProceduralGaussBlob(window, tw, th, backgroundcolorOffset);

% Aspect ratio width vs. height:
aspectratio = 1.0;

% Spatial constant
sc = 50;

%% Stimulus Sequence
% Flip before the loop
if ~isnan(DURATION_TRIGGER)
    % Send a trigger signal (98)
    sendtrigger(address, 98, DURATION_TRIGGER);
end

vbl = Screen('Flip', window);

%% Ramp up contrast
for n = 1:48 % rising side of taper
    lum = finTaper(n);
    Screen('DrawTexture', window, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [lum, sc, aspectratio, 0]);
    % Flip to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end

%% trigger
if ~isnan(DURATION_TRIGGER)
    % Send a trigger signal (99) after ramp
    sendtrigger(address, 55, DURATION_TRIGGER);
end

% Actual presentation interval:
for i = 1:numCycles 

    % if ~isnan(DURATION_TRIGGER)
    %     % Send a trigger signal (111)
    %     sendtrigger(address, 111, DURATION_TRIGGER);
    % end

    %% Contrast: -100
    contrast = -100;
    Screen('DrawTexture', window, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [contrast, sc, aspectratio, 0]);
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi); % Flip
    
    %% Contrast: 0
    contrast = 0;
    Screen('DrawTexture', window, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [contrast, sc, aspectratio, 0]);
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

    %% Contrast: 100
    contrast = 100;
    Screen('DrawTexture', window, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [contrast, sc, aspectratio, 0]);
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

    %% Contrast: 0
    contrast = 0;
    Screen('DrawTexture', window, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [contrast, sc, aspectratio, 0]);
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

end % loop 1:numCycles 

if ~isnan(DURATION_TRIGGER)
    % Send a trigger signal (211)
    sendtrigger(address, 211, DURATION_TRIGGER);
end

%% Ramp down contrast
for n = 49:96 % falling side of taper
    lum = finTaper(n);
    Screen('DrawTexture', window, blobtex, [], [], [], [], [], [], [], kPsychDontDoRotation, [lum, sc, aspectratio, 0]);
    % Flip to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end

if ~isnan(DURATION_TRIGGER)
    sendtrigger(address, 66, DURATION_TRIGGER);
end