function [pulsetrain] = pulseTrainF0(pitchMat, fs, hoptime, deg)

range = 40;
clip = range - 1;
% One f0 estimate per frame
% t = (0:100) * 0.01;
% f0 = ls10(100, 400, length(t));
f0 = vertcat(0, pitchMat(:,3), 0);
[numFrames, ~] = size(f0);
t = (0:numFrames-1) * hoptime;
f0 = f0';
% Interpolate to one estimate per sample
ti = min(t) : 1/fs : max(t);
f0i = interp1(t, f0, ti);

% Phasor spinning around at the appropriate rate
phase = cumsum(f0i * deg / fs);  % Maybe 4 * pi...

% Could make a pure sine wave like this
% sinEx = sin(phase);

% Actual excitation (soft-clipped sine wave)
pulsetrain = sigmoid(range * sin(phase) - clip);

% figure(1)
% plot(pulsetrain)
% waitforbuttonpress();
% Just for display purposes
% subplots(db(stft(ex, 1024, 1024, 256)))
end

function p = sigmoid(x)
    p = 1 ./ (1 + exp(-x));
end