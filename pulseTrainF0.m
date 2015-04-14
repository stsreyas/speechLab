function [pulsetrain, pitchMat] = pulseTrainF0(pitchMat, fs, hoptime, deg, type)

range = 400;
clip = range - 1;
% One f0 estimate per frame
% t = (0:100) * 0.01;
% f0 = ls10(100, 400, length(t));
pitchVec1 = vertcat(0, pitchMat(:,3), 0);
% try to interpolate when a single frame is out of place

pitchVec2 = removeJitter(pitchVec1);
pitchMat = pitchVec2(2:end-1);

pitchVec = spreadVec(pitchVec2, 4);

% figure(1)
% subplot(3,1,1)
% plot(pitchVec1)
% subplot(3,1,2)
% plot(pitchVec2)
% subplot(3,1,3)
% plot(pitchVec)
% waitforbuttonpress();

f0 = pitchVec;
[numFrames, ~] = size(f0);
t = (0:numFrames-1) * hoptime;
f0 = f0';
% Interpolate to one estimate per sample
ti = min(t) : 1/fs : max(t);
f0i = interp1(t, f0, ti);

% gauss1D = [1, 4, 6, 4, 1]'/16;
% f0iSmooth = imfilter(f0i, gauss1D, 'replicate');


% Phasor spinning around at the appropriate rate
phase = cumsum(f0i * deg / fs);  % Maybe 4 * pi...

% Could make a pure sine wave like this
% sinEx = sin(phase);
if strcmp(type, 'zc');
    pulsetrainzc = (diff(sin(phase)>0)>0);
    pulsetrain = pulsetrainzc;
% figure1 = figure;
% subplot(2,1,1)
% plot(sin(phase));
% subplot(2,1,2)
% plot(pulsetrainzc);
% waitforbuttonpress();
% close(figure1)

% Actual excitation (soft-clipped sine wave)
else
    pulsetrainsigmoid = sigmoid(range * sin(phase) - clip);
    pulsetrain = pulsetrainsigmoid;
end

% figure2 = figure;
% subplot(2, 1, 1)
% plot(pulsetrainsigmoid)
% subplot(2,1,2)
% plot(pulsetrainzc)
% waitforbuttonpress()
% close(figure2);

end

function p = sigmoid(x)
    p = 1 ./ (1 + exp(-x));
end

function output = removeJitter(input)

output = zeros(size(input));
[r, ~] = size(input);
startIdx = 3;
endIdx = r-2;

for i = startIdx:endIdx
   if(input(i) == 0)
       prevFrPitch = input(i-2);
       nextFrPitch = input(i+2);
       if((prevFrPitch > 0) && (nextFrPitch > 0))
            output(i) = (prevFrPitch + nextFrPitch) / 2;
       end
   else
       output(i) = input(i);
   end
end
end

function output = spreadVec(input, shift)

output = zeros(size(input));
[r, ~] = size(input);

startIdx = shift + 1;
endIdx = r - shift;

for i = startIdx:endIdx
    curFrPitch = input(i);
   if(curFrPitch ~= 0)
       prevFrPitch = input(i - shift);
       nextFrPitch = input(i + shift);
       if(prevFrPitch == 0)
            output((i - shift):i) = curFrPitch;
       elseif(nextFrPitch == 0)
            output(i:(i + shift)) = curFrPitch;
       else
           output(i) = input(i);
       end
%    else
%        output(i) = input(i);
   end
end
end