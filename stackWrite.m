function stackWrite(stack,varargin)
%STACKWRITE この関数の概要をここに記述
%   singleはuintにしてから

% singleかつラベルの対処
if isa(stack,'single') && max(stack(:))>1 || isa(stack,'double') && max(stack(:))>1
    nzStack=unique(stack);
    nzStack(nzStack==0)=[];
    minStack=min(nzStack(:));
    stack=stack-minStack+1;
    stack(stack<0)=0;
	stack=stack/max(stack(:));
    stack=im2uint16(stack);
end

if isa(stack,'logical')
    stack=mat2gray(stack);
end

if isempty(varargin)
    dirName='.\stack';
%     dirName='C:\Users\Azuma\Desktop\stack';
    fileName='img';
elseif length(varargin)==2
    dirName=varargin{1};
    fileName=varargin{2};
elseif length(varargin)==1
    dirName=varargin{1};
    fileName='img';
end
mkdir(dirName);

if ndims(stack)==4
    [r,c,z,t]=size(stack);
    for thisT=1:t
        for thisZ=1:z
            fileAddress=[dirName,'\',fileName,'_t',num2str(thisT,'%03u'),'_z',num2str(thisZ,'%03u'),'.tif'];
            imwrite(stack(:,:,thisZ,thisT),fileAddress,'tif');
        end
    end
elseif ndims(stack)==3
    [r,c,z]=size(stack);
    for thisZ=1:z
        fileAddress=[dirName,'\',fileName,'_z',num2str(thisZ,'%03u'),'.tif'];
        imwrite(stack(:,:,thisZ),fileAddress,'tif');
    end
elseif ismatrix(stack)
    fileAddress=[dirName,'\',fileName,'.tif'];
    imwrite(stack,fileAddress,'tif');
end
end

