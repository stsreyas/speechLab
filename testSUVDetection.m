pitchMat = textread('./data/440_16k/01/pitch/440c0201_noisy_airport.SAcC.pitch');
pitchMat2 = textread('./data/440_16k/01/pitch/440c0201_clean.SAcC.pitch');

voicedVec = (pitchMat(:,3) == 0)';
unvoicedVec = (pitchMat(:,4) < 0.4)';
silenceVec = horzcat(0, unvoicedVec .* voicedVec, 0);

[numFrames, ~] = size(silenceVec);

sStrVec = strfind(silenceVec, [0,1]);
sEndVec = strfind(silenceVec, [1,0]);

[~, nstarts] = size(sStrVec);
[~, nends] = size(sStrVec);

numV = min(nstarts, nends);
sCtr = 1;
eCtr = 1;

for i = 1:numV
    sStart = sStrVec(sCtr);
    sEnd = sEndVec(eCtr);
    if sEnd <= sStart
        sCtr = sCtr + 1;
        eCtr = eCtr + 1;
    else
        if (sEnd - sStart) < 20
            silenceVec(sStart:sEnd) = 0;
        end
        sCtr = sCtr + 1;
        eCtr = eCtr + 1;
    end
end

sStrVec = strfind(silenceVec, [0,1]);
sEndVec = strfind(silenceVec, [1,0]);

[~, nstarts] = size(sStrVec);
[~, nends] = size(sStrVec);

numV = min(nstarts, nends);
sCtr = 1;
eCtr = 1;

for i = 1:numV
    sStart = sStrVec(sCtr);
    sEnd = sEndVec(eCtr);
    if sEnd <= sStart
        sCtr = sCtr + 1;
        eCtr = eCtr + 1;
    else
        if i > 1
            prevSEnd = sEndVec(eCtr - 1);
            if(sStart - prevSEnd) < 20
                silenceVec(prevSEnd:sStart) = 1;
            end
        end
        sCtr = sCtr + 1;
        eCtr = eCtr + 1;
    end
end


figure(1)
subplot(2,1,1)
plot(pitchMat2(:,3));
subplot(2,1,2)
plot(silenceVec, 'b*-');