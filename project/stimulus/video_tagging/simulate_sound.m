% Parameters
audioLen = 17;          % total length (seconds)
nChirps = 15;            % number of chirps
chirp_dur = nChirps / audioLen - 0.2;        % duration of each chirp (seconds)
fs = 16000;             % sample rate
noise_level = 0.02;     % background noise amplitude

f_start = 400;          % chirp start freq (Hz)
f_end   = 4000;         % chirp end freq (Hz)
fadeDur = 0.04;         % fade in/out (seconds)
fadeSamp = round(fadeDur * fs);

% --- 1. Make background noise
audio = noise_level * randn(audioLen*fs, 1);

% --- 2. Insert chirps at regular intervals
chirp_gap = (audioLen - nChirps*chirp_dur)/(nChirps+1); % sec between chirps
chirp_onsets = round( (chirp_gap + (0:(nChirps-1)) * (chirp_dur + chirp_gap)) * fs ) + 1;
chirp_samps = round(chirp_dur * fs);
for k = 1:nChirps
    t = (0:chirp_samps-1)'/fs;
    w = chirp(t, f_start, chirp_dur, f_end, 'linear');   % w is a column
    w = 0.5 * w;                                         % scale

    fadeSamp_chirp = min(round(fadeDur*fs), length(w));
    fadein = linspace(0,1,fadeSamp_chirp)';              % also column
    fadeout = linspace(1,0,fadeSamp_chirp)';             % also column
    w(1:fadeSamp_chirp) = w(1:fadeSamp_chirp) .* fadein;
    w(end-fadeSamp_chirp+1:end) = w(end-fadeSamp_chirp+1:end) .* fadeout;

    idx1 = chirp_onsets(k);
    idx2 = min(idx1 + chirp_samps - 1, audioLen*fs);
    len = idx2 - idx1 + 1;
    audio(idx1:idx2) = audio(idx1:idx2) + w(1:len);
end

% --- 3. Normalize to avoid clipping
audio = audio / max(abs(audio)) * 0.98;

% --- 4. Save
outdir = fullfile(pwd, 'project', 'stimulus', 'audio_files');
if ~exist(outdir, 'dir'), mkdir(outdir); end

outfile = fullfile(outdir, ...
    sprintf('chirpNoise_%dchirps_%.2fsDur_%dHz-%dHz_%ds.wav', ...
        nChirps, chirp_dur, round(f_start), round(f_end), audioLen));

audiowrite(outfile, audio, fs);

% --- 5. Plot
t_full = (0:length(audio)-1)/fs;
figure; plot(t_full, audio, 'k');
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Chirp Audio (%d chirps in %ds, with noise)', nChirps, audioLen));
xlim([0 audioLen]); grid on;

% --- 6. Play
% sound(audio, fs)
fprintf('Audio file saved to: %s\n', outfile);
