function [filterMat] = filterbank(minFreq, maxFreq, numBins, fftSize, fs)

	totalSize = (fftSize / 2) + 1;
	filterMat = zeros(numBins, totalSize);
	freqCenters = getFreqCenters(minFreq, maxFreq, numBins, fftSize, fs);

	freqCenterBins = floor(freqCenters * fftSize / fs);

	for k = 1 : totalSize
		for m = 2:numBins+1
			if ((k >= freqCenterBins(m - 1)) && (k <= freqCenterBins(m)))
				numerator = 2 * (k - freqCenterBins(m - 1));
				denominator = (freqCenterBins(m + 1) - freqCenterBins(m - 1)) * (freqCenterBins(m) - freqCenterBins(m - 1));
				filterMat((m - 1),k) = numerator / denominator;
			
			elseif ((k >= freqCenterBins(m)) && (k <= freqCenterBins(m + 1)))
				numerator = 2 * (freqCenterBins(m + 1) - k);
				denominator = (freqCenterBins(m + 1) - freqCenterBins(m - 1)) * (freqCenterBins(m + 1) - freqCenterBins(m));
				filterMat((m - 1),k) = numerator / denominator;

			end
		end		
	end
end 

function [freqCenters] = getFreqCenters(minFreq, maxFreq, numBins, fftSize, fs)

	minB = 1125 * log( 1 + minFreq/700);
	maxB = 1125 * log( 1 + maxFreq/700);
	freqCenters = minFreq * ones((numBins + 2),1);

	for i = 2 : (numBins + 1)
		b = (minB + (i - 1) * ((maxB - minB)/(numBins + 1)));
		freqCenters(i) = 700 * (exp(b/1125) - 1);	
	end
	freqCenters(numBins + 2) = maxFreq;
end
