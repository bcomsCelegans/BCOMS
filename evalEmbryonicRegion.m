function [ overVal, vol, nucCoveredFlag ] = evalEmbryonicRegion( stack, thisNormImgSlice, thisNuc )
%EVALEMBREG この関数の概要をここに記述
%   詳細説明をここに記述

peri = bwperim(logical(stack), 8);

% 細胞膜として認識された領域のmembチャネルでの輝度値（高いほうが一致している）
if max(peri(:)) == 0
    overVal = 0;
else
    embRegOver=immultiply(thisNormImgSlice, peri);%評価するときはオリジナル画像
    overVal = mean(nonzeros(embRegOver));
end

% 体積
vol = sum(stack(:));

% 全細胞核が胚領域に含まれていれば１
[r,c,zNum] = size(thisNuc);
stack = reshape(stack, r, c, zNum);
nucOver1 = immultiply(thisNuc, ~stack);
if max(nucOver1(:)) == 1
    nucCoveredFlag = 0;
else
    nucCoveredFlag = 1;
end
