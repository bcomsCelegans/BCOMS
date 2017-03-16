function idx=idxN(stack)

id=find(stack);
sz=size(stack);

[out{1:ndims(stack)}] = ind2sub(sz,id);

idx = cell2mat(out);
