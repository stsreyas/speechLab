% Take the FFT of the window
% Create a mel filterbank for each window
% Apply filterbank
% Take the discrete cosine transform in the mel space
% Use the first 13 coefficients as the MFCC features

function [mfccVec] = mfcc(frame, fBank, numMFCC)

	[numBins, ~] = size(fBank);
	frameSamples = length(frame);
	fftSize = 2^nextpow2(frameSamples);
	x = fft(frame, fftSize);
	X = x(1 : (fftSize/2 + 1));
	XMagSq = ((real(X).^2) + (imag(X).^2));
	XMagSqMat = repmat(XMagSq', [numBins, 1]);
	S = zeros(numBins, 1);
	S = log(sum(XMagSqMat .* fBank, 2));
	c = discreteCosineTransform(S);
	mfccVec = c(1:numMFCC)';
end


function c = discreteCosineTransform(vec)

	c = zeros(size(vec));
	[numCoeff, ~] = size(vec);
	for i = 1 : numCoeff
		temp = 0;
		for j = 1 : numCoeff
			temp = temp + (vec(j) * cos(pi * i * (j -1/2) / numCoeff));
		end		
		c(i) = temp;
	end
end
