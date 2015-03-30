function [excitation, mfccPSpec] = generateExcitation(numSamples, type, step, varargin)

[sumpower, minfreq, maxfreq, nbands, bwidth, dcttype, fbtype, ...
    inPSpec, inMfccVec, amplitude, fs, nfft, pitch, hoptime] = ...
    process_options(varargin, 'sumpower', 0, 'minfreq', 50, 'maxfreq', 7000, ...
    'nbands', 26, 'bwidth', 1.0, 'dcttype', 1, 'fbtype', 'mel',...
    'inPSpec', [], 'inMfccVec', [], 'amp', 1, 'fs', 16000, 'nfft', 512, 'pitch', [], 'hoptime', 0.01);

excitation = [];
mfccPSpec = [];
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
    inMfccVec = lifter(inMfccVec, -22, 1);
    mfccSpec = cep2spec(inMfccVec, nbands, dcttype);
    mfccPSpec = invaudspec(mfccSpec, fs, nfft, fbtype, minfreq, maxfreq, sumpower, bwidth);
    excitationSpec = abs(inPSpec) ./ sqrt(mfccPSpec);
    excitationSpec(find(excitationSpec == Inf)) = inPSpec(find(excitationSpec == Inf));
    
%     figure1 = figure;
%     axes1 = axes('Parent', figure1);
%     hold(axes1, 'all');
%     subplot(3,1, 1)
%     imagesc(10*log10(inPSpec)); axis xy; colorbar    
%     subplot(3,1,2)
%     imagesc(10*log10(mfccPSpec)); axis xy; colorbar    
%     subplot(3,1,3)
%     imagesc(10*log10(excitationSpec)); axis xy; colorbar    
%     close(figure1);
    
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
    excitation(find(excitation == 0)) = noise;
end

end