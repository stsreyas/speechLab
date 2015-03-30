clc;
clear all;
addpath('E:\SLATE\Repos\rastamat\trunk\');
saveFiles = 1;

% noisyMFCC = load('./data/440c0201_noisy.mat');
% cleanMFCC = load('./data/440c0201_clean.mat');
% [cleanAud, ~] = audioread('./data/440c0201_clean.wav');
[cleanAud, ~] = audioread('./data/440_16k/01/440c0201_clean.wav');
[noisyAud, ~] = audioread('./data/440_16k/01/440c0201_noisy_airport.wav');


estMFCC = load('./data/440_16k/01/440c0201_it05_airport.mat');
pitchMat = textread('./data/440_16k/01/pitch/440c0201_clean.SAcC.pitch');
[numFrames2, ~] = size(pitchMat);

% Parameters being used in anBySyn
absOptions = { 'wintime', 0.020, 'hoptime', 0.010, 'sumpower', 1*0, 'preemph',  0.97*0,...
    'dither', 0, 'minfreq' ,50, 'maxfreq', 7000, 'bwidth', 1.0, 'modelorder', 0, 'nbands', 26,...
    'usecmp', 0, 'fbtype', 'htkmel', 'dcttype', 1, 'lifterexp', -22};

exc = 'suvpitchmix';
if strcmp(exc, 'residual'), excFlag = 1; end

fs = estMFCC.fs;
wintime = 0.020;
steptime = 0.010;
winpts = round(wintime * fs);
steppts = round(steptime * fs);
[~, numFrames1] = size(estMFCC.cep3);
numFrames = min(numFrames1, numFrames2);
numSamples = winpts + steppts*(numFrames- 1);

% Take mfcc vec
% return all the way to power spec
% divide the clean audio power spec with the mfcc ( both calculated on the
% clean audio and the estimated
% invert power spec to get two excitation models (ground truth &
% estimated)

moreOptions = {'amp', 10, 'fs', fs, 'nfft', 512, 'winpts', winpts,...
    'steppts', steppts, 'pitch', pitchMat};
absOptions = horzcat(absOptions, moreOptions);
ex = generateExcitation(numSamples, exc, steppts, absOptions{:});

% figure(1)
% plot(ex)
% waitforbuttonpress();

% excit = [];
excit = ex;

% [cep, ~, ~] = melfcc(aud, fs);
% [audResyn, ~, ~] = invmelfcc(cep, fs);
absOptions = horzcat(absOptions, {'excitation', excit});

[cleanAudPSpec, ~] = powspec(cleanAud, fs, wintime, steptime, 0);
[cleanAudMfcc, ~, ~] = melfcc(cleanAud, fs, absOptions{:});
[noisyAudPSpec, ~] = powspec(noisyAud, fs, wintime, steptime, 0);
[noisyAudMfcc, ~, ~] = melfcc(noisyAud, fs, absOptions{:});

estCep13 = estMFCC.cep3(1:13, 1:numFrames);

for i = 1:13
    
    maxClean = max(cleanAudMfcc(i, :));
    minClean = min(cleanAudMfcc(i, :));
    rangeAud = maxClean - minClean;
    
    maxNoisy = max(estCep13(i, :));
    minNoisy = min(estCep13(i, :));
    
    estCep13(i, :) = (estCep13(i, :) - minNoisy) / (maxNoisy- minNoisy);
    estCep13(i, :) = estCep13(i, :) * rangeAud + minClean;
    
end

[estResyn, estResynAspec, estResynPspec] = invmelfcc(estCep13, fs, 1, 0, absOptions{:});
[noisyResyn, noisyResynAspec, noisyResynPspec] = invmelfcc(noisyAudMfcc(:, 1:numFrames), fs, 1, 0, absOptions{:});

% sound(estResyn);

if saveFiles == 1
    set = '440c0201_it05_pc';
    fname = strcat(set, '_', exc, '.wav');
    audiowrite(fname, estResyn, fs);
    set = '440c0201_noisy_airport_pc';
    fname = strcat(set, '_', exc, '.wav');
    audiowrite(fname, noisyResyn, fs);
end

status = 3

figure(4)
subplot(2,1,1)
plot(cleanAud);
subplot(2,1,2)
plot(estResyn);

figure(5)
subplot(2,1,1)
plot(cleanAudMfcc(1, :));
subplot(2,1,2)
plot(estCep13(1, :));

figure(3)
subplot(3, 1, 1)
imagesc(10*log10(cleanAudPSpec)); axis xy; colorbar
subplot(3, 1, 2)
imagesc(10*log10(estResynAspec)); axis xy; colorbar
subplot(3, 1, 3)
imagesc(10*log10(estResynPspec)); axis xy; colorbar