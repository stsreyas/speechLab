function testPEASS
addpath('E:\SLATE\Repos\PEASS-Software-v2.0');
addpath('E:\SLATE\Repos\PEASS-Software-v2.0\gammatone');

testNum = 1;

% cleanFName = 'E:\SLATE\Repos\speechLab\trunk\data\440_16k\01\440c0201_clean.wav';
cleanFName = 'E:\SLATE\Repos\speechLab\trunk\440c0201_cleanResyn_residual.wav';
noisyFName = 'E:\SLATE\Repos\speechLab\trunk\440c0201_noisy_airport_residual.wav';
noisyDMFName = 'E:\SLATE\Repos\speechLab\trunk\440c0201_noisy_airportDM_DM.wav';
% estFName = 'E:\SLATE\Repos\speechLab\trunk\440c0201_it04_airport_residual.wav';
estFName = 'E:\SLATE\Repos\speechLab\trunk\440c0201_it04_residual2.wav';

if testNum == 1
    fnames = {cleanFName, estFName};
elseif testNum == 2
    fnames = {cleanFName, noisyFName};
elseif testNum == 3
    fnames = {cleanFName, noisyDMFName};
end

fnames = normalizeFileLength(fnames);

%%%%%%%%%%%%
% Set inputs
%%%%%%%%%%%%
originalFiles = fnames(1);
estimateFile = fnames{2};

%%%%%%%%%%%%%
% Set options
%%%%%%%%%%%%%
options.destDir = 'E:/SLATE/Repos/speechLab/trunk/peassRes/';
options.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

%%%%%%%%%%%%%%%%%%%%
% Call main function
%%%%%%%%%%%%%%%%%%%%
res = PEASS_ObjectiveMeasure(originalFiles,estimateFile,options);

%%%%%%%%%%%%%%%%%
% Display results
%%%%%%%%%%%%%%%%%

fprintf('************************\n');
fprintf('* INTERMEDIATE RESULTS *\n');
fprintf('************************\n');

fprintf('The decomposition has been generated and stored in:\n');
cellfun(@(s)fprintf(' - %s\n',s),res.decompositionFilenames);

fprintf('The ISR, SIR, SAR and SDR criteria computed with the new decomposition are:\n');
fprintf(' - SDR = %.1f dB\n - ISR = %.1f dB\n - SIR = %.1f dB\n - SAR = %.1f dB\n',...
    res.SDR,res.ISR,res.SIR,res.SAR);

fprintf('The audio quality (PEMO-Q) criteria computed with the new decomposition are:\n');
fprintf(' - qGlobal = %.3f\n - qTarget = %.3f\n - qInterf = %.3f\n - qArtif = %.3f\n',...
    res.qGlobal,res.qTarget,res.qInterf,res.qArtif);

fprintf('*************************\n');
fprintf('****  FINAL RESULTS  ****\n');
fprintf('*************************\n');
fprintf(' - Overall Perceptual Score: OPS = %.f/100\n',res.OPS)
fprintf(' - Target-related Perceptual Score: TPS = %.f/100\n',res.TPS)
fprintf(' - Interference-related Perceptual Score: IPS = %.f/100\n',res.IPS)
fprintf(' - Artifact-related Perceptual Score: APS = %.f/100\n',res.APS);
end

function fnames2 = normalizeFileLength(fnames)

[~, numFiles] = size(fnames);

minSize = 10^10;
for i = 1:numFiles
    [sig, ~] = audioread(fnames{i});
    [ns, ~] = size(sig);
    if ns < minSize
        minSize = ns;
    end
end

fnames2 = {};
for i = 1:numFiles
    [sig, fs] = audioread(fnames{i});
    fname2 = strcat(fnames{i}(1:end-4), '_sh.wav');
    audiowrite(fname2, sig(1:minSize), fs);
    fnames2 = horzcat(fnames2, fname2);
end

end
