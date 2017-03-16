function stack = tikan( stack, original, new )

[newStack,id] = ismember(stack(:),original(:)) ;
stack(newStack) = new(id(newStack)) ;