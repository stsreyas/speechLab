clc;
clear all;
addpath('E:\SLATE\Repos\rastamat\trunk\');
saveFiles = 1;

[cleanAud, fs] = audioread('./data/440_16k/01/440c0201_clean.wav');
[noisyAud, fs] = audioread('./data/440_16k/01/440c0201_noisy_airport.wav');

% Parameters being used in anBySyn
absOptions = { 'wintime', 0.020, 'hoptime', 0.010, 'sumpower', 1*0, 'preemph',  0.97*0,...
    'dither', 0, 'minfreq' ,50, 'maxfreq', 7000, 'bwidth', 1.0, 'modelorder', 0, 'nbands', 26,...
    'usecmp', 0, 'fbtype', 'mel', 'dcttype', 1, 'lifterexp', -22};

exc = 'residual';
pitchMat = [];
if strcmp(exc, 'residual'), excFlag = 1; end

wintime = 0.020;
steptime = 0.010;
winpts = round(wintime * fs);
steppts = round(steptime * fs);
nOverlap = winpts - steppts;
nfft = 512;

% [cleanAudPSpec, ~] = powspec(cleanAud, fs, wintime, steptime, 0);
hanningWindow = hanning(winpts)';
[cleanAudPSpec, ~] = specgram(cleanAud, nfft, fs, hanningWindow, nOverlap);
[cleanAudMfcc, ~, ~] = melfcc(cleanAud, fs, absOptions{:});
[noisyAudPSpec, ~] = specgram(noisyAud, nfft, fs, hanningWindow, nOverlap);
[noisyAudMfcc, ~, ~] = melfcc(noisyAud, fs, absOptions{:});

[~, numFrames] = size(cleanAudMfcc);
numSamples = winpts + steppts*(numFrames- 1);

moreOptions = {'amp', 10, 'fs', fs, 'nfft', 512, 'winpts', winpts,...
    'steppts', steppts, 'pitch', pitchMat, 'inPSpec', cleanAudPSpec, 'inMfccVec', cleanAudMfcc};
absOptions = horzcat(absOptions, moreOptions);

ex = generateExcitation(numSamples, exc, steppts, absOptions{:});

excit = ex;
absOptions = horzcat(absOptions, {'excitation', excit});

estMFCC = load('E:\SLATE\Repos\speechLab\trunk\data\440_16k\01\440c0201_it04_airport.mat');
% estMFCC13 = estMFCC.cep3(1:13, :);
estMFCC13 = estMFCC.state;

for i = 1:13
    
    maxClean = max(cleanAudMfcc(i, :));
    minClean = min(cleanAudMfcc(i, :));
    rangeClean = maxClean - minClean;
    
    maxEst01 = max(estMFCC13(i, :));
    minEst01 = min(estMFCC13(i, :));
    
    estMFCC13(i, :) = (estMFCC13(i, :) - minEst01) / (maxEst01- minEst01);
    estMFCC13(i, :) = estMFCC13(i, :) * rangeClean + minClean;
    
end

[cleanResyn, cleanResynASpec, cleanResynPSpec] = invmelfcc(cleanAudMfcc, fs, 1, 1, absOptions{:});
[noisyResyn, noisyResynASpec, noisyResynPSpec] = invmelfcc(noisyAudMfcc, fs, 1, 1, absOptions{:});
[estResyn, estResynASpec, estResynPSpec] = invmelfcc(estMFCC13, fs, 1, 1, absOptions{:});

if saveFiles == 1
    set = '440c0201_cleanResyn';
    fname = strcat(set, '_', exc, '.wav');
    audiowrite(fname, cleanResyn, fs);
    set = '440c0201_it04_airport';
    fname = strcat(set, '_', exc, '.wav');
    audiowrite(fname, estResyn, fs);
    set = '440c0201_noisy_airport';
    fname = strcat(set, '_', exc, '.wav');
    audiowrite(fname, noisyResyn, fs);
end

status = 3

figure(4)
subplot(4,1,1)
plot(cleanAud);
subplot(4,1,2)
plot(cleanResyn);
subplot(4,1,3)
plot(estResyn);
subplot(4,1,4)
plot(noisyResyn);

figure(5)
subplot(4,1,1)
imagesc(10*log10(abs(cleanAudPSpec).^2)); axis xy; colorbar
subplot(4,1,2)
imagesc(10*log10(abs(cleanResynPSpec).^2)); axis xy; colorbar
subplot(4,1,3)
imagesc(10*log10(abs(estResynPSpec).^2)); axis xy; colorbar
subplot(4,1,4)
imagesc(10*log10(abs(noisyResynPSpec).^2)); axis xy; colorbar