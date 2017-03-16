function [ overVal, vol, nucCoveredFlag ] = evalEmbryonicRegion( stack, thisNormImgSlice, thisNuc )
%EVALEMBREG ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q

peri = bwperim(logical(stack), 8);

% �זE���Ƃ��ĔF�����ꂽ�̈��memb�`���l���ł̋P�x�l�i�����ق�����v���Ă���j
if max(peri(:)) == 0
    overVal = 0;
else
    embRegOver=immultiply(thisNormImgSlice, peri);%�]������Ƃ��̓I���W�i���摜
    overVal = mean(nonzeros(embRegOver));
end

% �̐�
vol = sum(stack(:));

% �S�זE�j����̈�Ɋ܂܂�Ă���΂P
[r,c,zNum] = size(thisNuc);
stack = reshape(stack, r, c, zNum);
nucOver1 = immultiply(thisNuc, ~stack);
if max(nucOver1(:)) == 1
    nucCoveredFlag = 0;
else
    nucCoveredFlag = 1;
end
