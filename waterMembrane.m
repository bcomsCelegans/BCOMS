function h = waterMembrane( membImgDir, nucValDir, embRegDir, membSegDir, resXY, resZ, h )

% Temp
waterStackTempDir = [membSegDir, '\SegmentationTemp'];
mkdir(waterStackTempDir);
waterScoreTempDir = [membSegDir, '\ScoreTemp'];
mkdir(waterScoreTempDir);

% �œK��
waterStackDir = [membSegDir, '\Cell'];
mkdir(waterStackDir);
waterLabelDir = [membSegDir, '\CellLabel'];
mkdir(waterLabelDir);
waterMembDir = [membSegDir, '\Membrane'];
mkdir(waterMembDir);
matStackDir = [membSegDir, '\MatFile'];
mkdir(matStackDir);
% waterScoreDir = [membSegDir, '\Score'];
% mkdir(waterScoreDir);

% waitbar
waitbar(0.2, h);

%% �f�[�^����

% {

mkdir(membSegDir);

% Read
memb = oneStackLoad(membImgDir);

nuc = oneStackLoad(nucValDir);
nuc=double(nuc);

embReg=oneStackLoad(embRegDir);
embReg=double(logical(embReg));

% ����
ratio=round(resZ/resXY);
[~,~,zNum,tNum] = size(memb);


%% �O����
oriEmbReg = embReg;
embReg=imerode(embReg,ones(3, 3));
% embReg = arrayfun(@(x) imerode(embReg(:,:,:,x), ones(15-2*x)), 1:length(tList), 'UniformOutput', false);
% embReg = cat(4, embReg{:});
embReg(:,:,1,:) = embReg(:,:,1,:) * 0;
embReg(:,:,end,:) = embReg(:,:,end,:) * 0;

% �S�~����
embArea = sum(embReg(:)) / 10;
embReg = bwareaopen(embReg, round(embArea*0.1),26);

% �����ϊ�
dNuc = ratioDst( nuc, ratio, 3 );
dNuc(isinf(dNuc))=1;
dNuc=immultiply(double(dNuc),embReg);
dNuc=mat2gray(dNuc);

% memb��ave filt������܂ł̕��@��
% {
filtSize=[1 3 5];
% h=ones(filtSize)/prod(filtSize);
% membDenoise=imfilter(memb,h,'replicate');
%}

% waitbar
waitbar(0.3, h);

%% ������U����watershed�̎��s
alphas=0:0.005:0.04;
% alphas=0:0.1:0.4;
% alphas=0:4;

% {
% poolobj = gcp('nocreate');
% delete(poolobj);
% poolobj = parpool(2);
parfor si = 1:length(filtSize)
% for si = 1:length(filtSize)
    sz = filtSize(si);
    fs = [sz sz sz];
    hs=ones(fs)/prod(fs);
    membDenoise=imfilter(memb,hs,'replicate');
    for a=1:length(alphas)
    %     a = 2
        stack=membDenoise;
        alpha=alphas(a);

        stack=mat2gray(stack);

        stack=stack+alpha*dNuc;
        stack=mat2gray(stack);

        range = getrangefromclass(stack);

        % ������EmbReg��1st trial
        stack(~embReg)=range(2);

        % �זE�j�̈���ŏ��l�ɂ���
        stack=imimposemin(stack,logical(nuc),26);

        stack=watershed(stack,26);
        stack=immultiply(double(stack),embReg);

        % ID���j��ID�ɍ��킹��B�G�b�W�̈�Əd�Ȃ��Ă���j������̂łO�Ƃ̃y�A�͖�������
        conv = unique([stack(nuc>0) nuc(nuc>0)], 'rows');
        zeroLines = any(conv==0, 2);
        conv = conv(~zeroLines, :);
        stack=tikan(stack, conv(:,1), conv(:,2));

        membStack = logical(imdilate(logical(stack),ones(3,3,3)) - logical(stack));
        membOver=immultiply(membDenoise, membStack);%�_���p
%         membRevOver=immultiply(membDenoise, ~membStack);%�_���p

        % �P�x�l���O�łȂ��s�N�Z���̕��ϒl
        membOverMean = arrayfun(@(x) mean(nonzeros(membOver(:,:,:,x))), 1:size(membOver, 4));
%         memRevbOverMean = arrayfun(@(x) mean(nonzeros(membRevOver(:,:,:,x))), 1:size(membRevOver, 4));
        
        % s/n
%         membOverSN = membOverMean./memRevbOverMean;

        % �Ή�����זE�j��������Ȃ������זE�̐����e�^�C���|�C���g�Ōv�Z
        missedNum = arrayfun(@(x) length(setdiff(unique(nonzeros(nuc(:,:,:,x))), unique(nonzeros(stack(:,:,:,x))))), 1:tNum);

        % �e�זE�̐ς̍ő�ƍŏ��̔䗦
        uniColors = unique(nonzeros(stack));
        vol = arrayfun(@(x) squeeze(sum(sum(sum(stack==x, 1), 2), 3)), uniColors, 'UniformOutput', false);
        vol = cat(2, vol{:});
        vol = vol';%�s���זE��, �񂪎���
        vol(vol == 0) = NaN;

        volRatioAtTime = nanmin(vol) ./ nanmax(vol);


        score =  [(1:tNum)', repmat(sz, [length(membOverMean) 1]), repmat(alpha, [length(membOverMean) 1]), membOverMean', missedNum' repmat(volRatioAtTime, [length(membOverMean) 1])];
%         score =  [(1:tNum)', repmat(sz, [length(membOverSN) 1]), repmat(alpha, [length(membOverSN) 1]), membOverSN', missedNum' volRatioAtTime'];

        % �ۑ�
        filename = [waterStackTempDir, '\FiltSize', num2str(sz), '_Alpha', num2str(alpha), '.mat'];
        parsaveStack(filename, stack);
        filename = [waterScoreTempDir, '\FiltSize', num2str(sz), '_Alpha', num2str(alpha), '.mat'];
        parsaveScore(filename, score);
    end
end
%}

% waitbar
waitbar(0.8, h);

score = [];
for sz=filtSize
    for a=1:length(alphas)
        alpha=alphas(a);
        filename = [waterScoreTempDir, '\FiltSize', num2str(sz), '_Alpha', num2str(alpha), '.mat'];
        thisScore = oneStackLoad(filename);
%         thisScore = thisScore(4,:);
        fs = thisScore(1,2);
        alpha = thisScore(1,3);
        meanOver = mean(thisScore(:,4));
        numMissed = sum(thisScore(:,5));
        meanVolRatio = mean(thisScore(:,6));
        score = [score; fs, alpha, meanOver, numMissed, meanVolRatio];
    end
end

%% temp
% stack = [];
% for sz=filtSize
%     for a=1:length(alphas)
%         alpha=alphas(a);
%         filename = [waterStackTempDir, '\FiltSize', num2str(sz), '_Alpha', num2str(alpha), '.mat'];
%         thisStack = oneStackLoad(filename);
%         stack = cat(3, stack, thisStack(:,:,9,2));
% %         stack = cat(3, stack, thisStack(:,:,20,4));
%     end
% end
% stack = randomizeId(stack);
% % visualize4Dsc(stack);
% memb = [];
% for sz=filtSize
%     for a=1:length(alphas)
%         alpha=alphas(a);
%         filename = [waterStackTempDir, '\FiltSize', num2str(sz), '_Alpha', num2str(alpha), '.mat'];
%         thisStack = oneStackLoad(filename);
%         thisStack = thisStack(:,:,9,2);
%         thisMemb = logical(immultiply(oriEmbReg(:,:,9,2), ~logical(thisStack)));
%         memb = cat(3, memb, thisMemb);
%     end
% end
% stack = stack + (max(stack(:) + 1)) * logical(memb);

%% �S�������̓K�p
% {
% missed �זE��
numMissedThre = 0;
score = score(score(:,4) == numMissedThre,:);

%}

%% �ړI�֐��̕]��
% �œK�����v�Z
score = sortrows(score, 3);
minSz = score(end,1);
minAlpha = score(end,2);

% Load
filename = [waterStackTempDir, '\FiltSize', num2str(minSz), '_Alpha', num2str(minAlpha), '.mat'];
stack = oneStackLoad(filename);

rmdir(waterStackTempDir, 's');
rmdir(waterScoreTempDir, 's');

stack = bwconncomp(stack, 26);
stack = labelmatrix(stack);

% waitbar
waitbar(0.9, h);

%% �ۑ�
% stack = randomizeId(stack);
stackWrite(stack, waterLabelDir);
stackWrite(logical(stack), waterStackDir);

filename = [matStackDir, '\cell.mat'];
parsaveStack(filename, stack)

% stackWrite(memb, waterMembDir);

% saveName=[waterScoreDir,'\score.mat'];
% parsaveScore(saveName, score);

% % with Membrane
% memb = logical(immultiply(oriEmbReg, ~logical(stack)));
% stackWithMemb = randomizeId(stack);
% memb = stackWithMemb + (max(stackWithMemb(:) + 1)) * logical(memb);
% memb(:,:,1,:) = memb(:,:,1,:) * 0;
% memb(:,:,end,:) = memb(:,:,end,:) * 0;
% visualize4Dsc(memb);

% Only Membrane
memb = logical(immultiply(oriEmbReg, ~logical(stack)));
memb(:,:,1,:) = bwperim(memb(:,:,1,:), 8);
memb(:,:,end,:) = bwperim(memb(:,:,end,:), 8);
memb(:,:,1,:) = memb(:,:,1,:) * 0;
memb(:,:,end,:) = memb(:,:,end,:) * 0;
stackWrite(memb, waterMembDir);

