clc;
clear all;
addpath('E:\SLATE\Repos\rastamat\trunk\');
saveFiles = 1;

[cleanAud, fs] = audioread('./data/440_16k/01/440c0201_clean.wav');
[noisyAud, fs] = audioread('./data/440_16k/01/440c0201_noisy_airport.wav');
estMFCC = load('E:\SLATE\Repos\speechLab\trunk\data\440_16k\01\440c0201_it04_airport.mat');
estMask = estMFCC.mask;

% Parameters being used in anBySyn
% absOptions = { 'wintime', 0.020, 'hoptime', 0.010, 'sumpower', 1*0, 'preemph',  0.97*0,...
%     'dither', 0, 'minfreq' ,50, 'maxfreq', 7000, 'bwidth', 1.0, 'modelorder', 0, 'nbands', 26,...
%     'usecmp', 0, 'fbtype', 'mel', 'dcttype', 1, 'lifterexp', -22};

sumpower = 1*0; 
preemph = 0.97*0;
dither = 0;
minfreq = 50;
maxfreq = 7000;
bwidth = 1.0;
modelorder = 0;
nbands = 26;
usecmp = 0;
fbtype = 'mel';
dcttype = 1;
lifterexp = -22;

wintime = 0.020;
steptime = 0.010;
winpts = round(wintime * fs);
steppts = round(steptime * fs);
nOverlap = winpts - steppts;
nfft = 512;

hanningWindow = hanning(winpts)';
[cleanAudPSpec, ~] = specgram(cleanAud, nfft, fs, hanningWindow, nOverlap);
[noisyAudPSpec, ~] = specgram(noisyAud, nfft, fs, hanningWindow, nOverlap);

estMaskPSpec = invaudspec(estMask, fs, nfft, fbtype, minfreq, maxfreq, sumpower, bwidth);

cleanAudPSpecMasked = sqrt(estMaskPSpec) .* cleanAudPSpec;
noisyAudPSpecMasked = sqrt(estMaskPSpec) .* noisyAudPSpec;
% [cleanResynMasked] = invspecgram(cleanAudPSpecMasked, nfft, fs, winpts, (winpts - steppts));
% [noisyResynMasked] = invspecgram(noisyAudPSpecMasked, nfft, fs, winpts, (winpts - steppts));
cleanResynMasked = ispecgram(cleanAudPSpecMasked, nfft, fs, winpts, (winpts - steppts));
noisyResynMasked = ispecgram(noisyAudPSpecMasked, nfft, fs, winpts, (winpts - steppts));


if saveFiles == 1
    set = '440c0201_cleanResyn';
    fname = strcat(set, '_DM.wav');
    audiowrite(fname, cleanResynMasked, fs);
    set = '440c0201_noisy_airport';
    fname = strcat(set, '_DM.wav');
    audiowrite(fname, noisyResynMasked, fs);
end

figure(4)
subplot(2,1,1)
plot(cleanResynMasked);
subplot(2,1,2)
plot(noisyResynMasked);

figure(5)
subplot(5,1,1)
imagesc(10*log10(abs(estMaskPSpec).^2)); axis xy; colorbar
subplot(5,1,2)
imagesc(10*log10(abs(cleanAudPSpec).^2)); axis xy; colorbar
subplot(5,1,3)
imagesc(10*log10(abs(cleanAudPSpecMasked).^2)); axis xy; colorbar
subplot(5,1,4)
imagesc(10*log10(abs(noisyAudPSpec).^2)); axis xy; colorbar
subplot(5,1,5)
imagesc(10*log10(abs(noisyAudPSpecMasked).^2)); axis xy; colorbar