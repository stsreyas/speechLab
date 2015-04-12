function [excitation, mfccPSpec] = generateExcitation(numSamples, type, step, varargin)

[sumpower, minfreq, maxfreq, nbands, bwidth, dcttype, fbtype, ...
    inPSpec, inMfccVec, amplitude, fs, nfft, pitch, hoptime, winpts, steppts] = ...
    process_options(varargin, 'sumpower', 0, 'minfreq', 50, 'maxfreq', 7000, ...
    'nbands', 26, 'bwidth', 1.0, 'dcttype', 1, 'fbtype', 'mel',...
    'inPSpec', [], 'inMfccVec', [], 'amp', 1, 'fs', 16000, 'nfft', 512,...
    'pitch', [], 'hoptime', 0.01, 'winpts', 400, 'steppts', 160);

excitation = [];
mfccPSpec = [];
probThresh = 0.6;
% if nargin < 10
% if nargin < 9, dcttype = 1;
% if nargin < 8, nOverlap = 160; end
% if nargin < 7, winpts = 400; end
% if nargin < 6, fs = 16000; end
% if nargin < 5, nfft = 512; end

if strcmp(type, 'pulsetrain')
    excitation = zeros([numSamples, 1]);
    excitation(1:step:numSamples) = amplitude;
    
elseif strcmp(type, 'gaussianpulse')
    excitation = zeros([numSamples, 1]);
    x = (-ceil(step/2):1:floor(step/2));
    norm = 10000*normpdf(x,0,1);
    numReps = floor(numSamples/step);
    excitation = repmat(norm', [numReps,1]);
    
elseif strcmp(type, 'whitenoise')
    excitation = zeros([numSamples, 1]);
    excitation = randn(numSamples,1);
    
elseif strcmp(type, 'residual')
    if (size(inPSpec(:,1)) <= 1), error('input signal power spectrum not present'); end
    if (size(inMfccVec(:,1)) <= 1), error('mfcc power spectrum not present'); end
    % lifter is applicable here since cep3 is the liftered mfcc values
    inMfccVec = lifter(inMfccVec, -22, 1);
    mfccSpec = cep2spec(inMfccVec, nbands, dcttype);
    mfccPSpec = invaudspec(mfccSpec, fs, nfft, fbtype, minfreq, maxfreq, sumpower, bwidth);
%     excitationSpec = sqrt(inPSpec) ./ sqrt(mfccPSpec);
    excitationSpec = (inPSpec) ./ sqrt(mfccPSpec);
%     excitationSpec(excitationSpec == Inf) = inPSpec(excitationSpec == Inf);
    absExcitationSpec = abs(excitationSpec);
    excitationSpec(absExcitationSpec == Inf) = 0;
    excitationSpec(225:end, :) = 0;
    excitationSpec(1:2, :) = 0;
    
%     tempSpec = excitationSpec .* sqrt(mfccPSpec);
%     [tempSig] = invspecgram(tempSpec, nfft, fs, winpts, (winpts - steppts));
%     
%     saveFiles = 1;
%     if saveFiles == 1
%         set = '440c0201_tempResyn';
%         fname = strcat(set, '_', type, '.wav');
%         audiowrite(fname, tempSig, fs);
%     end
% 
%     figure1 = figure;
%     axes1 = axes('Parent', figure1);
%     hold(axes1, 'all');
%     subplot(3,1,1)
%     imagesc(10*log10(abs(inPSpec).^2)); axis xy; colorbar    
%     subplot(3,1,2)
%     imagesc(10*log10(mfccPSpec)); axis xy; colorbar    
%     subplot(3,1,3)
%     imagesc(10*log10(abs(tempSpec).^2)); axis xy; colorbar    
%     waitforbuttonpress();
    
%     nOverlap = winpts - steppts;
%     excitation = ispecgram(excitationSpec, nfft, fs, winpts, nOverlap);
    excitation = excitationSpec;
    
elseif strcmp(type, 'pitchpulsetrain')
    excitation = pulseTrainF0(pitch, fs, hoptime, (2 * pi));

elseif strcmp(type, 'pitchnoisemix')
    pulseEx = pulseTrainF0(pitch, fs, hoptime, (2 * pi));
    [ns, ~] = size(pulseEx);
    noise = randn(ns,1);
    excitation = pulseEx;
    excitation(pulseEx == 0) = noise;
    
elseif strcmp(type, 'suv')
    
    suvPVec = vertcat(0, pitch(:,4));
    suvVec = suvPVec > probThresh;
    suvMat = repmat(suvVec', [(winpts-steppts), 1]);
    [r, c] = size(suvMat);
    suvSamples = reshape(suvMat, [1, r*c]);
    pulseEx = zeros([numSamples, 1]);
    pulseEx(1:step:numSamples) = amplitude;
    noiseEx = randn(numSamples,1);
    excitation = zeros(size(suvSamples));
    excitation(suvSamples == 0) = noiseEx(suvSamples == 0);
    excitation(suvSamples == 1) = pulseEx(suvSamples == 1);
    
elseif strcmp(type, 'suvpitchmix')
    
    suvPVec = vertcat(0, pitch(:,4));
    suvVec = suvPVec > probThresh;
    suvMat = repmat(suvVec', [(winpts-steppts), 1]);
    [r, c] = size(suvMat);
    suvSamples = reshape(suvMat, [1, r*c]);

    [pulseEx, ~] = pulseTrainF0(pitch, fs, hoptime, (2 * pi), 'zc' );
    noiseEx = randn(numSamples,1);
    excitation = zeros(size(suvSamples));
    excitation(suvSamples == 0) = noiseEx(suvSamples == 0);
    excitation(suvSamples == 1) = pulseEx(suvSamples == 1);

elseif strcmp(type, 'suvpitchmix2')
    
    scaleNoise = 1;
    scalePitch = 10;
    [pulseEx, modPitchVec] = pulseTrainF0(pitch, fs, hoptime, (2 * pi), 'zc');
    noiseEx = rand(numSamples,1);

    pitch(:,3) = modPitchVec;
    suvPVec = vertcat(0, pitch(:,3));
    suvVec = suvPVec == 0;
    suvMat = repmat(suvVec', [(winpts-steppts), 1]);
    [r, c] = size(suvMat);
    suvSamples = reshape(suvMat, [1, r*c]);
    suvSamples = abs(suvSamples - 1);
%     figure(1)
%     plot(suvSamples, 'b*');
%     waitforbuttonpress();

    % temporary cross fade 
    excitation = zeros(size(suvSamples));
    excitation(suvSamples == 0) = scaleNoise * noiseEx(suvSamples == 0);% - scaleNoise/2;
    excitation(suvSamples == 1) = scalePitch * pulseEx(suvSamples == 1);% - scalePitch/2;
    
end

end