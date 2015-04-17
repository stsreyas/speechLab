[cleanAud, fs] = audioread('E:\SLATE\Repos\speechLab\trunk\data\440_16k\01\440c0201_clean.wav');

% [cleanAud, fs] = audioread('E:/Dropbox/results/residualResyn/440c0201_cleanResyn_residual.wav');
% [estAud, ~] = audioread('E:/Dropbox/results/residualResyn/440c0201_it04_airport_residual.wav');
[estAud, ~] = audioread('E:\SLATE\Repos\speechLab\trunk\440c0201_it04_residual2.wav');

% [noisyAud, ~] = audioread('E:/Dropbox/results/residualResyn/440c0201_noisy_airport_residual.wav');
[noisyAud, ~] = audioread('E:\SLATE\Repos\speechLab\trunk\440c0201_noisy_airportDM_DM.wav');
% [cleanAud, fs] = audioread('E:/SLATE/Repos/speechLab/trunk/pitchRes/440c0201_clean_zc_suvpitchmix2_zc_1lpf_1hpf_decay_5.wav');
% [estAud, ~] = audioread('E:/SLATE/Repos/speechLab/trunk/pitchRes/440c0201_it04_pn_zc_suvpitchmix2_zc_1lpf_1hpf_decay_5.wav');
% [noisyAud, ~] = audioread('E:/SLATE/Repos/speechLab/trunk/pitchRes/440c0201_noisy_airport_zc_suvpitchmix2_zc_1lpf_1hpf_decay_5.wav');

% sound(cleanAud, fs);
% pause(5);
% sound(noisyAud, fs);
% pause(5);
% sound(estAud, fs);
% pause(5);

[numSamples1, ~] = size(estAud);
[numSamples2, ~] = size(cleanAud);

numSamples = min(numSamples1, numSamples2);

dEst = stoi(cleanAud(1:numSamples, :), estAud(1:numSamples, :), fs);
dNoisy = stoi(cleanAud(1:numSamples, :), noisyAud(1:numSamples, :), fs);

status = 1