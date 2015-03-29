clc;
clear all;
addpath('E:\SLATE\Repos\rastamat\trunk\');

saveFiles = 1;
% Parameters being used in anBySyn
excVec = {'whitenoise', 'pulsetrain', 'residual',  'gaussianpulse'};
ampVec = [0.5, 5, 50];
winVec = [0.020, 0.025];
lifterExpVec = [0.6, -22];
lifterFlagVec = [0, 1];
dirVec = {'\01\', '\05\'};

goAgain = 1;
[~, numExc] = size(excVec);
[~, numAmp] = size(ampVec);
[~, numWinSize] = size(winVec);
[~, numLifterExp] = size(lifterExpVec);
[~, numLifterFlags] = size(lifterFlagVec);
[~, numDirs] = size(dirVec);
dirRoot = 'E:\SLATE\Repos\speechLab\trunk\data\440_16k';

for dirIter = 1:numDirs
    dName = strcat(dirRoot, dirVec(dirIter));
    directoryName = dName{:};
    namingConventionClean = 'clean';
    namingConventionNoisy = 'noisy';
    namingConventionEstimate = 'it';
    
    wildcardClean = strcat(directoryName, '/*', namingConventionClean, '*');
    fileListingClean = dir(wildcardClean);
    wildcardNoisy = strcat(directoryName, '/*', namingConventionNoisy, '*');
    fileListingNoisy = dir(wildcardNoisy);
    wildcardEstimate = strcat(directoryName, '/*', namingConventionEstimate, '*');
    fileListingEstimate = dir(wildcardEstimate);
    
    [numClean, ~] = size(fileListingClean);
    [numNoisy, ~] = size(fileListingNoisy);
    [numEstimate, ~] = size(fileListingEstimate);
    
    % need to make this more general
    [clean, ~] = audioread(strcat(directoryName, '/', fileListingClean(1).name));
    [noisy, ~] = audioread(strcat(directoryName, '/', fileListingNoisy(1).name));
    estimate01 = load(strcat(directoryName, '/', fileListingEstimate(1).name));
    estimate02 = load(strcat(directoryName, '/', fileListingEstimate(2).name));
    fs = estimate01.fs;
    [~, numFrames] = size(estimate01.cep3);
    steptime = 0.010;
    excFlag = 0;
    
    for excIter = 1:numExc-1
        exc = excVec(excIter);
        if strcmp(exc, 'residual'), excFlag = 1; end
        for ampIter = 2:numAmp
            amp = ampVec(ampIter);
            for winIter = 2:numWinSize
                wintime = winVec(winIter);
                for lifterExpIter = 2:numLifterExp
                    lifterexp = lifterExpVec(lifterExpIter);
                    for lifterFlagIter = 2:numLifterFlags
                        lifterflag = lifterFlagVec(lifterFlagIter);
                        winpts = round(wintime * fs);
                        steppts = round(steptime * fs);
                        numSamples = winpts + steppts*(numFrames- 1);
                        
                        mainOpts = { 'wintime', wintime, 'hoptime', 0.010, 'sumpower', 0, 'preemph',  0,...
                            'dither', 0, 'minfreq', 50, 'maxfreq', 7000, 'bwidth', 1.0, 'modelorder', 0, 'nbands', 26,...
                            'usecmp', 0, 'fbtype', 'mel', 'dcttype', 1, 'lifterexp', lifterexp,...
                            'amp', amp, 'fs', fs, 'nfft', 512, 'winpts', winpts, 'steppts', steppts};
                        
                        % Take mfcc vec
                        % return all the way to power spec
                        % divide the clean audio power spec with the mfcc ( both calculated on the
                        % clean audio and the estimated
                        % invert power spec to get two excitation models (ground truth &
                        % estimated)
                        
                        [cleanPSpec, ~] = powspec(clean, fs, wintime, steptime, 0);
                        [noisyPSpec, ~] = powspec(noisy, fs, wintime, steptime, 0);
                        
                        [cleanMfcc, ~, ~] = melfcc(clean, fs, mainOpts{:});
                        [noisyMfcc, ~, ~] = melfcc(noisy, fs, mainOpts{:});
                        
                        cleanOptions = {'inPSpec', cleanPSpec, 'inMfccVec', cleanMfcc};
                        %                         noisyOptions = {'inPSpec', noisyPSpec, 'inMfccVec', noisyMfcc};
                        
                        cleanOptionsFull = horzcat(mainOpts, cleanOptions);
                        %                         noisyOptionsFull = horzcat(mainOpts, noisyOptions);
                        
                        [cleanEx, cleanMfccPSpec] = generateExcitation(numSamples, exc, steppts, cleanOptionsFull{:});
                        %                         noisyEx = generateExcitation(numSamples, exc, steppts, noisyOptionsFull{:});
                        
                        cleanOptionsFull = horzcat(cleanOptionsFull, {'excitation', cleanEx});
                        %                         noisyOptionsFull = horzcat(noisyOptionsFull, {'excitation', noisyEx});
                        
                        [cleanResyn, cleanResynASpec, cleanResynPSpec] = ...
                            invmelfcc(cleanMfcc(1:13,:), fs, lifterflag, excFlag, cleanOptionsFull{:});
                        [noisyResyn, noisyResynASpec, noisyResynPSpec] = ...
                            invmelfcc(noisyMfcc(1:13,:), fs, lifterflag, excFlag, cleanOptionsFull{:});
                        [~, numFramesClean] = size(cleanMfcc);
                        
                        estimate01Cep13 = estimate01.cep3(1:13, 1:numFramesClean);
                        estimate02Cep13 = estimate02.cep3(1:13, 1:numFramesClean);
                        
                        for i = 1:13
                            
                            maxClean = max(cleanMfcc(i, :));
                            minClean = min(cleanMfcc(i, :));
                            rangeClean = maxClean - minClean;
                            
                            maxEst01 = max(estimate01Cep13(i, :));
                            minEst01 = min(estimate01Cep13(i, :));
                            
                            maxEst02 = max(estimate02Cep13(i, :));
                            minEst02 = min(estimate02Cep13(i, :));
                            
                            estimate01Cep13(i, :) = (estimate01Cep13(i, :) - minEst01) / (maxEst01- minEst01);
                            estimate01Cep13(i, :) = estimate01Cep13(i, :) * rangeClean + minClean;
                            
                            estimate02Cep13(i, :) = (estimate02Cep13(i, :) - minEst02) / (maxEst02- minEst02);
                            estimate02Cep13(i, :) = estimate02Cep13(i, :) * rangeClean + minClean;
                            
                        end
                        
                        [est01Resyn, est01ResynASpec, est01ResynPSpec] = ...
                            invmelfcc(estimate01Cep13, fs, lifterflag, excFlag, cleanOptionsFull{:});
                        [est02Resyn, est02ResynASpec, est02ResynPSpec] = ...
                            invmelfcc(estimate02Cep13, fs, lifterflag, excFlag, cleanOptionsFull{:});
                        
                        if saveFiles == 1
                            folderName = strcat(directoryName, '\results');
                            mkdir(folderName);
                            cleanSet = fileListingClean(1).name(1:end-4);
                            cleanname = sprintf('%s\\%s_%s_amp%.1f_wt%.2f_le%.2f_lf%d.wav',...
                                folderName, cleanSet, exc{:}, amp, wintime, lifterexp, lifterflag);
                            audiowrite(cleanname, cleanResyn, fs);
                            
                            noisySet = fileListingNoisy(1).name(1:end-4);
                            noisyname = sprintf('%s\\%s_%s_amp%.1f_wt%.2f_le%.2f_lf%d.wav',...
                                folderName, noisySet, exc{:}, amp, wintime, lifterexp, lifterflag);
                            audiowrite(noisyname, noisyResyn, fs);
                            
                            est01Set = fileListingEstimate(1).name(1:end-4);
                            est01Name = sprintf('%s\\%s_%s_amp%.1f_wt%.2f_le%.2f_lf%d.wav',...
                                folderName, est01Set, exc{:}, amp, wintime, lifterexp, lifterflag);
                            audiowrite(est01Name, est01Resyn, fs);
                            
                            est02Set = fileListingEstimate(2).name(1:end-4);
                            est02Name = sprintf('%s\\%s_%s_amp%.1f_wt%.2f_le%.2f_lf%d.wav',...
                                folderName, est02Set, exc{:}, amp, wintime, lifterexp, lifterflag);
                            audiowrite(est02Name, est02Resyn, fs);
                        end
                        
                        %                         status = 3
                        
                        plotName1 = sprintf('%s\\resynSignals_%s_amp%.1f_wt%.2f_le%.2f_lf%d.jpg',...
                            folderName, exc{:}, amp, wintime, lifterexp, lifterflag);
                        
                        figure1 = figure;
                        axes1 = axes('Parent', figure1);
                        hold(axes1, 'all');
                        subplot(4,1,1)
                        plot(cleanResyn);
                        subplot(4,1,2)
                        plot(noisyResyn);
                        subplot(4,1,3)
                        plot(est01Resyn);
                        subplot(4,1,4)
                        plot(est02Resyn);
                        saveas(figure1, plotName1);
                        close(figure1);
                        
                        plotName2 = sprintf('%s\\audSpec_%s_amp%.1f_wt%.2f_le%.2f_lf%d.jpg',...
                            folderName, exc{:}, amp, wintime, lifterexp, lifterflag);
                        
                        figure2 = figure;
                        axes2 = axes('Parent', figure2);
                        hold(axes2, 'all');
                        subplot(4, 1, 1)
                        imagesc(10*log10(cleanResynASpec)); axis xy; colorbar
                        subplot(4, 1, 2)
                        imagesc(10*log10(noisyResynASpec)); axis xy; colorbar
                        subplot(4, 1, 3)
                        imagesc(10*log10(est01ResynASpec)); axis xy; colorbar
                        subplot(4, 1, 4)
                        imagesc(10*log10(est02ResynASpec)); axis xy; colorbar
                        saveas(figure2, plotName2);
                        close(figure2);
                        
                        plotName3 = sprintf('%s\\powerSpec_%s_amp%.1f_wt%.2f_le%.2f_lf%d.jpg',...
                            folderName, exc{:}, amp, wintime, lifterexp, lifterflag);
                        
                        figure3 = figure;
                        axes3 = axes('Parent', figure3);
                        hold(axes3, 'all');
                        subplot(4, 1, 1)
                        imagesc(10*log10(cleanResynPSpec)); axis xy; colorbar
                        subplot(4, 1, 2)
                        imagesc(10*log10(noisyResynPSpec)); axis xy; colorbar
                        subplot(4, 1, 3)
                        imagesc(10*log10(est01ResynPSpec)); axis xy; colorbar
                        subplot(4, 1, 4)
                        imagesc(10*log10(est02ResynPSpec)); axis xy; colorbar
                        saveas(figure3, plotName3);
                        close(figure3);
                        
                        if strcmp(exc, 'residual')
                            plotName4 = sprintf('%s\\idealExc_%s_amp%.1f_wt%.2f_le%.2f_lf%d.jpg',...
                                folderName, exc{:}, amp, wintime, lifterexp, lifterflag);
                            figure4 = figure;
                            axes4 = axes('Parent', figure4);
                            hold(axes4, 'all');
                            subplot(3,1, 1)
                            imagesc(10*log10(cleanPSpec)); axis xy; colorbar
                            subplot(3,1,2)
                            imagesc(10*log10(cleanMfccPSpec)); axis xy; colorbar
                            subplot(3,1,3)
                            imagesc(10*log10(cleanEx)); axis xy; colorbar
                            saveas(figure4, plotName4);
                            close(figure4);
                        end
                    end
                end
            end
        end
    end
end