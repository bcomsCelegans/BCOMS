function embryonicRegionOriginal(membImg, nucValDir, embRegDir, volRatioThresh)

mkdir(embRegDir);

embRegEvalTempDir=[embRegDir, '\EvalTemp'];
mkdir(embRegEvalTempDir);
embRegParamTempDir=[embRegDir, '\ParamTemp'];
mkdir(embRegParamTempDir);
embRegStackTempDir=[embRegDir, '\StackTemp'];
mkdir(embRegStackTempDir);
embRegStackDir=[embRegDir, '\Stack'];
mkdir(embRegStackDir);

% ç◊ñEñå
memb=oneStackLoad(membImg);

% ç◊ñEäj
nuc=oneStackLoad(nucValDir);
nuc = logical(nuc);

% nucÇ™ë∂ç›ÇµÇƒÇ¢ÇÈéûä‘ÇæÇØéÊÇËèoÇ∑
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

% ScoreåvéZÇÃÇΩÇﬂÇÃê≥ãKâª
membNorm = reshape(memb, [r, c, zNum*tNum]);
means = arrayfun(@(x) mean(reshape(membNorm(:,:,x), [r*c, 1])), 1:size(membNorm,3));
means = means / mean(means);
membNorm = arrayfun(@(x) membNorm(:,:,x) / means(x), 1:size(membNorm,3), 'UniformOutput', false);
membNorm = cat(3, membNorm{:});
membNorm = reshape(membNorm, [r, c, zNum, tNum]);
membNorm = membNorm / max(membNorm(:));

% ê≥ãKâªZ
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

% timeÇï™â
try
    p = parpool;
catch
    p = gcp;
end
colNumT = p.NumWorkers;
setNumFactor = 1;
colNumT = colNumT * setNumFactor;
rowNumT = ceil(tNum / colNumT);
% èáî‘ÇÉâÉìÉ_ÉÄÇ…ï¿Ç◊ë÷Ç¶ÇÈ
tListMod = randperm(tNum, tNum);
% ïsë´ï™ÇÕ0Çí«â¡Ç∑ÇÈ
addNum = colNumT*rowNumT - tNum;
if addNum > 0
    tListMod = [tListMod zeros(1,addNum)];
end
tListMod = reshape(tListMod, [rowNumT colNumT]);

param.sigma=[0 1.5 3];
param.smooth=0.5;
param.contBiasFactorA=[0.005: 0.005 0.03];
param.contBiasFactorB=1:3:5;
param.repeatFactor=50:50:300;

% param.sigma=[1.5];
% param.smooth=0.5;
% param.contBiasFactorA=[0.3];
% param.contBiasFactorB=[3];
% param.repeatFactor=[1000];

eval = []; eval.vol = []; eval.nucRejected = []; eval.objective = [];
% {
for rn=1:rowNumT
    thisT = tListMod(rn,:);
    param = structfun(@(x) unique(x), param, 'UniformOutput', false);
    % {
    parfor i = 1:length(thisT)
%     for i = 1:length(thisT)
        t = thisT(i);
        if t==0
            continue
        end
        
        thisMemb = memb(:,:,:,t);
        thisNuc = nuc(:,:,:,t);
        thisIniReg = iniReg(:,:,:,t);
        thisMembNorm = membNorm(:,:,:,t);
        thisIni = reshape(thisIniReg, [r, c*zNum]);
        thisNormImgSlice = reshape(thisMembNorm, [r, c*zNum]);
        tempScore = [];
        
        % èâä˙âª
        evalTemp = []; evalTemp.vol = []; evalTemp.nucRejected = []; evalTemp.objective = [];
        paramTemp = []; paramTemp.sigma = []; paramTemp.smooth = []; paramTemp.contBiasFactorA = []; paramTemp.contBiasFactorB = []; paramTemp.repeatFactor = [];
        
        for sigma=param.sigma
            if sigma~=0
                gausImg=imgaussfilt(thisMemb, sigma);
            else
                gausImg=thisMemb;
            end
            for sm = param.smooth
                for ca = param.contBiasFactorA
                    for cb = param.contBiasFactorB
                        contBias = ca * vars(t).^cb;
                        for rf = param.repeatFactor
                            thisImg = reshape(gausImg, [r, c*zNum]);
                            membReg=activecontour(thisImg, thisIni,rf,'Chan-Vese' ,'SmoothFactor',sm, 'ContractionBias',contBias);

                            % eval
                            [ overVal, vol, nucCoveredFlag ] = evalEmbryonicRegion( membReg, thisNormImgSlice, thisNuc );

                            % constraintsÇÃï]âø
                            evalTemp.vol = [evalTemp.vol vol];
                            evalTemp.nucRejected = [evalTemp.nucRejected ~nucCoveredFlag];

                            % objective functionÇÃï]âø
                            evalTemp.objective = [evalTemp.objective overVal];

                            % segmentation result
                            paramTemp.sigma = [paramTemp.sigma sigma];
                            paramTemp.smooth = [paramTemp.smooth sm];
                            paramTemp.contBiasFactorA = [paramTemp.contBiasFactorA ca];
                            paramTemp.contBiasFactorB = [paramTemp.contBiasFactorB cb];
                            paramTemp.repeatFactor = [paramTemp.repeatFactor rf];
                        end
                    end
                end
            end
        end
        filename = [embRegEvalTempDir, '\T', num2str(t), '.mat'];
        parsaveStack(filename, evalTemp);
        filename = [embRegParamTempDir, '\T', num2str(t), '.mat'];
        parsaveStack(filename, paramTemp);
    end
    %}
    % evalÇÇ‹Ç∆ÇﬂÇÈ
    for t = thisT
%         t
        if t==0
            continue
        end
        filename = [embRegEvalTempDir, '\T', num2str(t), '.mat'];
        evalTemp = oneStackLoad(filename);
        filename = [embRegParamTempDir, '\T', num2str(t), '.mat'];
        param = oneStackLoad(filename);

        % Nucelar enclosure
        eval.nucRejected(t,:) = evalTemp.nucRejected;

        % volume ratio
        eval.vol(t,:) = evalTemp.vol;

        % objective
        eval.objective(t,:) = evalTemp.objective;
    end

    % constraintsÇÃï]âø
    nucRejected = logical(sum(eval.nucRejected));
    
    volMin = min(eval.vol);
    volMax = max(eval.vol);
    volRatio = volMin ./ volMax;
    volRejected = volRatio < volRatioThresh;
    
    rejected = logical(nucRejected  + volRejected);
    
    % temp
    rejected = rejected*0;
    
    % paramÇÃçXêV
    param = structfun(@(x) x(:, ~rejected), param, 'UniformOutput', false);
    eval = structfun(@(x) x(:, ~rejected), eval, 'UniformOutput', false);
    
    % äeÉpÉâÉÅÅ[É^ÅAÇ∑Ç◊ÇƒÇÃílÇ≈rejectÇ≥ÇÍÇΩèÍçá
    emptyFlag = structfun(@(x) isempty(x), param);
    if any(emptyFlag)
        errordlg('No segmentation was acquired in the parameter space');
    end
end
%}

% ç≈ìKâÇÃíTçı
meanObj = mean(eval.objective, 1);
[~, maxCol] = max(meanObj);
optParam = structfun(@(x) x(maxCol), param, 'UniformOutput', false);

% ç≈ìKâÇ≈çƒåvéZ
for rn=1:rowNumT
    thisT = tListMod(rn,:);
    parfor i = 1:length(thisT)
%     for i = 1:length(thisT)
        t = thisT(i);
        if t==0
            continue
        end
        
        thisMemb = memb(:,:,:,t);
        thisIniReg = iniReg(:,:,:,t);
        thisIni = reshape(thisIniReg, [r, c*zNum]);

        for sigma=optParam.sigma
            if sigma~=0
                gausImg=imgaussfilt(thisMemb, sigma);
            else
                gausImg=thisMemb;
            end
            for sm = optParam.smooth
                for ca = optParam.contBiasFactorA
                    for cb = optParam.contBiasFactorB
                        contBias = ca * vars(t).^cb;
                        for rf = optParam.repeatFactor
                            thisImg = reshape(gausImg, [r, c*zNum]);
                            membRegOpt=activecontour(thisImg, thisIni,rf,'Chan-Vese' ,'SmoothFactor',sm, 'ContractionBias',contBias);
                        end
                    end
                end
            end
        end
        membRegOpt = reshape(membRegOpt, [r, c, zNum]);
        filename = [embRegStackTempDir, '\T', num2str(t), '.mat'];
        parsaveStack(filename, membRegOpt);
    end
    embReg = zeros(size(memb));
    for t = thisT
        if t==0
            continue
        end
        filename = [embRegStackTempDir, '\T', num2str(t), '.mat'];
        stackTemp = oneStackLoad(filename);
        embReg(:,:,:,t) = stackTemp;
    end
end


% ï€ë∂
filename = [embRegStackDir, '\embrayonicRegion.mat'];
parsaveStack(filename, embReg);

