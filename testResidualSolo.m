clc;
clear all;
addpath('E:\SLATE\Repos\rastamat\trunk\');
saveFiles = 1;

[cleanAud, fs] = audioread('./data/440_16k/01/440c0201_clean.wav');
sumpower = 1*0; 
preemph = 0.97*0;
dither = 0;
minfreq = 50;
maxfreq = 7000;
bwidth = 1.0;
modelorder = 0;
nbands = 26;
usecmp = 0;
fbtype = 'htkmel';
dcttype = 1;
lifterexp = -22;
% Parameters being used in anBySyn

exc = 'residual';
pitchMat = [];
if strcmp(exc, 'residual'), excFlag = 1; end

wintime = 0.020;
steptime = 0.010;
winpts = round(wintime * fs);
steppts = round(steptime * fs);
nOverlap = winpts - steppts;
nfft = 512;

absOptions = { 'wintime', wintime, 'hoptime', steptime, 'sumpower', sumpower, 'preemph',  preemph,...
    'dither', dither, 'minfreq' ,minfreq, 'maxfreq', maxfreq, 'bwidth', bwidth, 'modelorder', modelorder, 'nbands', nbands,...
    'usecmp', usecmp, 'fbtype', fbtype, 'dcttype', dcttype, 'lifterexp', lifterexp};

% [cleanAudPSpec, ~] = powspec(cleanAud, fs, wintime, steptime, 0);
hanningWindow = hanning(winpts)';
% [cleanAudPSpec, ~] = specgram(cleanAud*32768, nfft, fs, hanningWindow, nOverlap);
[cleanAudPSpec, ~] = specgram(cleanAud, nfft, fs, hanningWindow, nOverlap);
[cleanAudMfcc, ~, ~] = melfcc(cleanAud, fs, absOptions{:});

cleanAudMfcc = lifter(cleanAudMfcc, -22, 1);
mfccSpec = cep2spec(cleanAudMfcc, nbands, dcttype);
mfccPSpec = invaudspec(mfccSpec, fs, nfft, fbtype, minfreq, maxfreq, sumpower, bwidth);
%     excitationSpec = sqrt(inPSpec) ./ sqrt(mfccPSpec);
excitationSpec = (cleanAudPSpec) ./ sqrt(mfccPSpec);
%     excitationSpec(excitationSpec == Inf) = inPSpec(excitationSpec == Inf);
absExcitationSpec = abs(excitationSpec);
excitationSpec(absExcitationSpec == Inf) = 0;
excitationSpec(225:end, :) = 0;
excitationSpec(1:2, :) = 0;

cleanResynPSpec = excitationSpec .* sqrt(mfccPSpec);
[cleanResyn] = invspecgram(cleanResynPSpec, nfft, fs, winpts, (winpts - steppts));

if saveFiles == 1
    set = '440c0201_cleanResyn';
    fname = strcat(set, '_', exc, '.wav');
    audiowrite(fname, cleanResyn, fs);
end

status = 3

figure(4)
subplot(2,1,1)
plot(cleanAud);
subplot(2,1,2)
plot(cleanResyn);
figure(5)
subplot(2,1,1)
imagesc(10*log10(abs(cleanAudPSpec).^2)); axis xy; colorbar
subplot(2,1,2)
imagesc(10*log10(abs(cleanResynPSpec).^2)); axis xy; colorbar