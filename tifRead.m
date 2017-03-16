function [ output_args ] = tifRead( memraneFilename, zNum, tNum, imageDir )
%TIFREAD この関数の概要をここに記述
%   詳細説明をここに記述

info = imfinfo(memraneFilename);
num_images = numel(info);
r = info(1).Height;
c = info(1).Width;
imageStack = zeros(r, c, num_images);
for k = 1:num_images
    imageStack(:,:,k) = imread(memraneFilename, k);
    % ... Do something with image A ...
end
imageStack = reshape(imageStack, [r,c,zNum, tNum]);
mkdir(imageDir);
imageFilename = [imageDir, '\stack.mat'];
parsaveStack(imageFilename, imageStack);

