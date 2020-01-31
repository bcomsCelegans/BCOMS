function h = embryonicRegion(membImg, nucValDir, embRegDir, volRatioThresh, h)

mkdir(embRegDir);

embRegStackTempDir=[embRegDir, '\StackTemp'];
mkdir(embRegStackTempDir);
embRegStackDir=[embRegDir, '\Stack'];
mkdir(embRegStackDir);

% Membrane
memb=oneStackLoad(membImg);

% Nucleus
nuc=oneStackLoad(nucValDir);
nuc = logical(nuc);

% nuc existing times
nucPix = idxN(nuc);
if ndims(nuc) >= 4
    tList = min(nucPix(:,4)):max(nucPix(:,4));
else
    tList = 1;
end

% size
[r,c,zNum,~]=size(memb);
tNum = max(tList);

memb=memb(:,:,:,tList);
nuc=nuc(:,:,:,tList);

% initial contour
for i = 1:10
    sumMemb = sum(memb, 4);
    sumMembDenoise = imgaussfilt(sumMemb, 6);
    meanInt = mean(sumMembDenoise(:));
    mask = sumMembDenoise > meanInt - (meanInt * 0.1*(i-1));

    mask = reshape(mask, [r, c*zNum]);
    imgReshape = reshape(sumMembDenoise, [r, c*zNum]);

    smooth=0.05;
    contBias=0.01;
    repeatTime = 100;

    mask=activecontour(imgReshape, mask,repeatTime,'Chan-Vese' ,'SmoothFactor',smooth, 'ContractionBias',contBias);
    mask = reshape(mask, [r, c, zNum]);
    mask = imdilate(mask, ones(7,7));

    % initial contour evaluation
    outNuc = immultiply(~mask, sum(nuc, 4));
    if max(outNuc(:)) == 0
        break
    end
end

iniReg = repmat(mask, [1 1 1 tNum]);

if isnan(volRatioThresh)
    filename = [embRegStackDir, '\embrayonicRegion.mat'];
    parsaveStack(filename, iniReg);
    return
end

% Normalize for score calculation
membNorm = reshape(memb, [r, c, zNum*tNum]);
means = arrayfun(@(x) mean(reshape(membNorm(:,:,x), [r*c, 1])), 1:size(membNorm,3));
means = means / mean(means);
membNorm = arrayfun(@(x) membNorm(:,:,x) / means(x), 1:size(membNorm,3), 'UniformOutput', false);
membNorm = cat(3, membNorm{:});
membNorm = reshape(membNorm, [r, c, zNum, tNum]);
membNorm = membNorm / max(membNorm(:));

% Normalize Z
means = arrayfun(@(x) mean(reshape(memb(:,:,x,:), [r*c*tNum, 1])), 1:size(memb,3));
means = means / mean(means);
means(means==0) = 0.0001;
memb = arrayfun(@(x) memb(:,:,x,:) / means(x), 1:size(memb,3), 'UniformOutput', false);
memb = cat(3, memb{:});
memb = memb / max(memb(:));

vars = arrayfun(@(x) var(reshape(memb(:,:,:,x), [r*c*zNum, 1])), 1:size(memb,4));
vars = vars / mean(vars);
if isnan(vars)
    vars = 1;
end

% {

try
    p = parpool;
catch
    p = gcp;
end
colNumT = p.NumWorkers;
setNumFactor = 1;
colNumT = colNumT * setNumFactor;
rowNumT = ceil(tNum / colNumT);
% Randemize
tListMod = randperm(tNum, tNum);
% Add 0 for short
addNum = colNumT*rowNumT - tNum;
if addNum > 0
    tListMod = [tListMod zeros(1,addNum)];
end
tListMod = reshape(tListMod, [rowNumT colNumT]);


% waitbar
waitbar(0.2, h);

% Parameters
smooth=0.5;
contBiasFactorA = [0.03 0.05:0.05:0.2];%default
contBiasFactorB=[4 6 8 10 12];%default
repeatFactor=100;%default
% repeatFactor=50;%50
erdSz = [1 3 5];
% repeatFactor=50;%50
% Score
score = {};
paramLength = length(contBiasFactorA) * length(contBiasFactorB) * length(repeatFactor);
% {

% parfor t=1:tNum
for t=1:tNum
    thisMemb = memb(:,:,:,t);
    % Gaussian filter
    thisMemb=imgaussfilt3(thisMemb, 1);
    tempScore = zeros(paramLength, 8);
    repeatFactorUpdated = repeatFactor;
    thisMembNorm = membNorm(:,:,:,t);
    thisNucErd = nuc(:,:,:,t);
    i = 0;
    for ca=contBiasFactorA
        for cb=contBiasFactorB
            contBias = ca * vars(t).^cb;
            for rf = repeatFactorUpdated
                membReg = zeros(r,c,zNum);
                parfor z = 1:zNum
%                 for z = 1:zNum
                    thisMemb2D = thisMemb(:,:,z);
                    thisIniReg2D = iniReg(:,:,z,t);
                    thisMean = mean(thisMemb2D(:));
                    thisStd = std(thisMemb2D(:));
                    thisCC = thisStd / thisMean;
                    if isnan(thisCC)
                        thisCC = 1;
                    end
                    thisSmooth = smooth / thisCC^2;
                    thisContBias = contBias / thisCC^2;
                    thisRf = round(rf * thisCC);
                    % Level set
                    membReg2D=activecontour(thisMemb2D, thisIniReg2D,thisRf,'Chan-Vese' ,'SmoothFactor',thisSmooth, 'ContractionBias',thisContBias);
                    membReg(:,:,z) = membReg2D;
                end
                for er = erdSz
                    i = i + 1;
                    membReg = imerode(membReg, ones(er, er));
                    % optimization
                    peri = bwperim(logical(membReg), 8);
                    if max(peri(:)) == 0
                        overVal = 0;
                    else
                        embRegOver=immultiply(thisMembNorm, peri);%Evaluation on the original image
                        overVal = mean(nonzeros(embRegOver));
                    end
                    % vol
                    vol = sum(membReg(:));
                    % 1 if all nuclei emclosed
                    membReg = reshape(membReg, r, c, zNum);
                    nucOver = immultiply(thisNucErd, ~membReg);
                    if max(nucOver(:)) == 1
                        nucCoveredFlag = 0;
                    else
                        nucCoveredFlag = 1;
                    end

                    % Score
                    tempScore(i, :) = [t ca cb rf er overVal vol nucCoveredFlag];
                    % save
                    filename = [embRegStackTempDir, '\T', num2str(t), 'CA', num2str(ca), 'CB', num2str(cb), 'RF', num2str(rf), 'ER', num2str(er), '.mat'];
                    parsaveStack(filename, membReg);
                end
            end
        end
        score{t} = tempScore;
    end
end


% waitbar
waitbar(0.8, h);

%%Find optimal parameters
% {
% Average through the time
scoreStack = cat(3, score{:});
meanScore = mean(scoreStack, 3);
minScore = min(scoreStack, [], 3);
maxScore = max(scoreStack, [], 3);
ratioScore = minScore ./ maxScore;

% {
% Constraints
% Volume consistency
volRatioCol = 7;
% volRatioThresh = 0.93;
meanScore = meanScore(ratioScore(:,volRatioCol) >= volRatioThresh,:);

objCol = 6;
meanScore = sortrows(meanScore, objCol, 'descend' );

% Finish if no satisaction
if isempty(meanScore)
    h = msgbox('No segmentation result satisfied the volume ratio constraint');
    return
end

%% ç≈Recalculation on the optimal parameters
% Parameters
caOpt=meanScore(1,2);
cbOpt=meanScore(1,3);
rfOpt=meanScore(1,4);
optEr=meanScore(1,5);

embOpt = zeros(r,c,zNum, tNum);
for t = 1:tNum
    filename = [embRegStackTempDir, '\T', num2str(t), 'CA', num2str(caOpt), 'CB', num2str(cbOpt), 'RF', num2str(rfOpt), 'ER', num2str(optEr), '.mat'];
    embOpt(:,:,:,t) = oneStackLoad(filename);
end


% waitbar
waitbar(0.9, h);

% save
filename = [embRegStackDir, '\embrayonicRegion.mat'];
parsaveStack(filename, embOpt);

rmdir(embRegStackTempDir, 's')

