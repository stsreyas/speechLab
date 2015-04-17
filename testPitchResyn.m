clc;
clear all;
addpath('E:\SLATE\Repos\rastamat\trunk\');
saveFiles = 1;

folderName = 'E:\SLATE\Repos\speechLab\trunk\pitchRes\';

[cleanAud, ~] = audioread('./data/440_16k/01/440c0201_clean.wav');
[noisyAud, ~] = audioread('./data/440_16k/01/440c0201_noisy_airport.wav');
estMFCC = load('./data/440_16k/01/440c0201_it04_airport.mat');
% pitchMat = textread('./data/440_16k/01/pitch/440c0201_noisy_airport.SAcC.pitch');
pitchMat = textread('./data/440_16k/01/pitch/440c0201_clean.SAcC.pitch');
[numFrames2, ~] = size(pitchMat);

hpfFlagVec = [0,1];
lpfFlagVec = [0,1];
pulseTypeVec = {'zc', 'other'};
fadeTypeVec = {'step', 'decay'};
fadeFramesVec = 1:5;
fadeSamplesVec = 160:800;

% Parameters being used in anBySyn
absOptions = { 'wintime', 0.020, 'hoptime', 0.010, 'sumpower', 1*0, 'preemph',  0.97*0,...
    'dither', 0, 'minfreq', 50, 'maxfreq', 7000, 'bwidth', 1.0, 'modelorder', 0, 'nbands', 26,...
    'usecmp', 0, 'fbtype', 'htkmel', 'dcttype', 1, 'lifterexp', -22};

exc = 'suvpitchmix3';
if strcmp(exc, 'residual'), excFlag = 1; end

fs = estMFCC.fs;
wintime = 0.020;
steptime = 0.010;
nfft = 512;
winpts = round(wintime * fs);
steppts = round(steptime * fs);
nOverlap = winpts - steppts;
[~, numFrames1] = size(estMFCC.cep3);
numFrames = min(numFrames1, numFrames2);
numSamples = winpts + steppts*(numFrames- 1);

% Take mfcc vec
% return all the way to power spec
% divide the clean audio power spec with the mfcc ( both calculated on the
% clean audio and the estimated
% invert power spec to get two excitation models (ground truth &
% estimated)

moreOptions = {'amp', 10, 'fs', fs, 'nfft', nfft, 'winpts', winpts,...
    'steppts', steppts, 'pitch', pitchMat};
absOptions = horzcat(absOptions, moreOptions);

lpFilt = designfilt('lowpassiir','FilterOrder', 4, ...
    'PassbandFrequency', 6000,'PassbandRipple', .5, ...
    'SampleRate', fs);

hpFilt = designfilt('highpassiir','FilterOrder', 4, ...
    'PassbandFrequency' , 400, 'PassbandRipple', .5, ...
    'SampleRate', fs);

pulseType = pulseTypeVec{1};
for hpfIter = 1:2
    hpfFlag = hpfFlagVec(hpfIter);
    for lpfIter = 1:2
        lpfFlag = lpfFlagVec(lpfIter);
        for fadeTypeIter = 2:2
            fadeType = fadeTypeVec{fadeTypeIter};
            for fadeWidthIter = 1:5
                fadeFrames = fadeFramesVec(fadeWidthIter);
                fadeSamples = fadeSamplesVec(fadeWidthIter);
                pitchOptions = {'pulseType', pulseType, 'fadeFrames', fadeFrames,...
                    'fadeSamples', fadeSamples, 'fadeType', fadeType};
                absOptions = horzcat(absOptions, pitchOptions);
                ex = generateExcitation(numSamples, exc, steppts, absOptions{:});
                
                if lpfFlag
                    ex = filter(lpFilt, ex);
                end
                if hpfFlag
                    ex = filter(hpFilt, ex);
                end
                hanningWindow = hanning(winpts)';
                [exPSpec, ~] = specgram(ex, nfft, fs, hanningWindow, nOverlap);
                
                plotName1 = sprintf('%s\\excitation_%s_%s_%dlpf_%dhpf_%s_%d.jpg',...
                    folderName, exc, pulseType, lpfFlag, hpfFlag, fadeType, fadeFrames);
                figure1 = figure;
                subplot(2,1,1)
                plot(ex)
                subplot(2,1,2)
                imagesc(10*log10(abs(exPSpec).^2)); axis xy; colorbar
                saveas(figure1, plotName1);
                close(figure1);
                
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
                [cleanResyn, cleanResynAspec, cleanResynPspec] = invmelfcc(cleanAudMfcc(:, 1:numFrames), fs, 1, 0, absOptions{:});
                [noisyResyn, noisyResynAspec, noisyResynPspec] = invmelfcc(noisyAudMfcc(:, 1:numFrames), fs, 1, 0, absOptions{:});
                
                % sound(estResyn);
                
                if saveFiles == 1
                    set1 = '440c0201_it04_pn_zc';
                    estName = sprintf('%s\\%s_%s_%s_%dlpf_%dhpf_%s_%d.wav',...
                        folderName,set1, exc, pulseType, lpfFlag, hpfFlag, fadeType, fadeFrames);
                    audiowrite(estName, estResyn, fs);
                    set2 = '440c0201_noisy_airport_zc';
                    noisyName = sprintf('%s\\%s_%s_%s_%dlpf_%dhpf_%s_%d.wav',...
                        folderName, set2, exc, pulseType, lpfFlag, hpfFlag, fadeType, fadeFrames);
                    audiowrite(noisyName, noisyResyn, fs);
                    set3 = '440c0201_clean_zc';
                    cleanName = sprintf('%s\\%s_%s_%s_%dlpf_%dhpf_%s_%d.wav',...
                        folderName, set3, exc, pulseType, lpfFlag, hpfFlag, fadeType, fadeFrames);
                    audiowrite(cleanName, cleanResyn, fs);
                end
                
                plotName2 = sprintf('%s\\resynSignals_%s_%s_%dlpf_%dhpf_%s_%d.jpg',...
                    folderName, exc, pulseType, lpfFlag, hpfFlag, fadeType, fadeFrames);
                figure2 = figure;
                subplot(2,1,1)
                plot(cleanAud);
                subplot(2,1,2)
                plot(estResyn);
                saveas(figure2, plotName2);
                close(figure2);
                % figure(5)
                % subplot(2,1,1)
                % plot(cleanAudMfcc(1, :));
                % subplot(2,1,2)
                % plot(estCep13(1, :));
                
                plotName3 = sprintf('%s\\resynPSpec_%s_%s_%dlpf_%dhpf_%s_%d.jpg',...
                    folderName, exc, pulseType, lpfFlag, hpfFlag, fadeType, fadeFrames);
                figure3 = figure;
                subplot(2, 1, 1)
                imagesc(10*log10(cleanAudPSpec)); axis xy; colorbar
                % subplot(3, 1, 2)
                % imagesc(10*log10(estResynAspec)); axis xy; colorbar
                subplot(2, 1, 2)
                imagesc(10*log10(abs(estResynPspec).^2)); axis xy; colorbar
                saveas(figure3, plotName3);
                close(figure3);
            end
        end
    end
end