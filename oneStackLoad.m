function stack = oneStackLoad( preDir, varargin )
%ONEFILEOPEN この関数の概要をここに記述
%   詳細説明をここに記述

if ~isdir(preDir)
    stack=load(preDir);
    try
        stack=stack.stack;
    catch
        try
            stack=stack.data;
        catch
            stack=stack.score;
        end
    end
    return
end

if isempty(varargin)
    i = 1;
else
    i = varargin{1};
end

D=dir(preDir);
preFilenameArr=[];
for k=1:length(D)
    if strcmp(D(k).name,'.') || strcmp(D(k).name,'..')
        continue;
    end
    preFilename=[D(k).name];
    preFilenameArr=char(preFilenameArr,preFilename);
end
if ~isempty(preFilenameArr)
    preFilenameArr(1,:)=[];%Remove first blank
end
    
% for i=1:size(preFilenameArr, 1)
% for i=1
filename=[preDir, filesep, preFilenameArr(i,:)];
stack=load(filename);
try
    stack=stack.stack;
catch
    try
        stack=stack.data;
    catch
        stack=stack.score;
    end
end

% end
