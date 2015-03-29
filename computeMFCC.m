% Input the signal
% form frames
% take hamming weight
% call mfcc for each frame

function mfccMat = computeMFCC(x, fs, varargin)

% Other arguments are
%  winTime - the size in ms of each frame
%  stepTime - the jump in ms between each frame
%  applyHamming - 1 if hamming window is to be applier
%  numMFCC - number of mfcc's to be saved
%  numBins - number of bins in the filterbank

nargin

if nargin < 7, numBins = 24; end
if nargin < 6, numMFCC = 13; end
if nargin < 5, applyHamming = 1; end
if nargin < 4, stepTime = 0.010; end
if nargin < 3, winTime = 0.025; end
if nargin < 2, fs = 16000; end
if nargin < 1, error('input signal not present'); end

windowSamples = winTime * fs;
stepSamples = stepTime * fs;
totalSamples = length(x);
hammingWindow = hamming(windowSamples);
fftSize = 2^nextpow2(windowSamples);

minFreq = 0;
maxFreq = fs / 2;

fBank = filterbank(minFreq, maxFreq, numBins, fftSize, fs);

mfccMat = [];
for i = 1 : stepSamples : (totalSamples - windowSamples) 

	frame = x(i : (i + windowSamples - 1));
	if applyHamming == 1	
		winFrame = frame .* hammingWindow;
	end
	mfccVec = mfcc(winFrame, fBank, numMFCC);
	mfccMat = vertcat(mfccMat, mfccVec);
end

end

