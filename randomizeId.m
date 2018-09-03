function stack = randomizeId( stack, varargin )
%RANDOMIZEID stack��0�ȊO�̐��l�������_���ɕύX����
% varargin(1):0�ȊO�ɏ��O���鐔�l�̃��X�g
% varargin(2):�ύX��̐��l�̍ő�l

% if length(varargin)==2
%     rmv=varargin{1};
%     levelUp=varargin{2};
% elseif length(varargin)==1
%     rmv=varargin{1};
% else
%     rmv=[];
%     levelUp=[];
% end
% 
% if ~isempty(rmv)
%     
% end


% if ~isempty(varargin)
%     levelUp=varargin{1};
% else
%     levelUp=[];
% end

if length(varargin)==2
    levelUp=varargin{1};
    maxId=varargin{2};
elseif length(varargin)==1
    levelUp=varargin{1};
    maxId=[];
else
    levelUp=[];
    maxId=[];
end

uniStack=nonzeros(unique(stack));
if ~isempty(maxId)
    rnd=randperm(maxId-levelUp);
    rnd = repmat(rnd, [1, ceil(length(uniStack) / (maxId-levelUp))]);
    rnd = rnd(1:length(uniStack));
else
    rnd=randperm(length(uniStack));
end

if ~isempty(levelUp)
    rnd=rnd+levelUp;
end
stack=replace(stack,uniStack,rnd);

end

