function stack = ratioDst( stack, ratio, dim )
%RATIODST この関数の概要をここに記述
%   詳細説明をここに記述

if dim==4
    if ratio(1)==1 && ratio(2)==1
        stack=bwdist(stack);
    elseif ratio(1)>=1 && ratio(2)>=1
        [x,y,z,t]=size(stack);
        ratioStack=zeros(x,y,ratio(1)*z,ratio(2)*t);
        for thisT=1:t
            for thisZ=1:z
                ratioStack(:,:,ratio(1)*thisZ,ratio(2)*thisT)=stack(:,:,thisZ,thisT);
            end
        end
        ratioStack=bwdist(ratioStack);
        stack=zeros(x,y,z,t);
        for thisT=1:t
            for thisZ=1:z
                stack(:,:,thisZ,thisT)=ratioStack(:,:,ratio(1)*thisZ,ratio(2)*thisT);
            end
        end
    elseif ratio(1)<1 || ratio(2)<1
        ratioAll=[1,1,ratio(1),ratio(2)];
        if ratioAll(3)<1
            ratioAll=ratioAll/ratioAll(3);
        end
        if ratioAll(4)<1
            ratioAll=ratioAll/ratioAll(4);
        end
        
%         ratioStack=zeros(x*ratioAll(1),y*ratioAll(2),z*ratioAll(3),t*ratioAll(4));
        bigStack=imresize(stack,ratioAll(1),'nearest');
        [x,y,z,t]=size(bigStack);
        ratioStack=zeros(x,y,ratioAll(1)*z,ratioAll(2)*t);
        for thisT=1:t
            for thisZ=1:z
                ratioStack(:,:,ratioAll(1)*thisZ,ratioAll(2)*thisT)=bigStack(:,:,thisZ,thisT);
            end
        end
        ratioStack=bwdist(ratioStack);
        bigStack=zeros(x,y,z,t);
        for thisT=1:t
            for thisZ=1:z
                bigStack(:,:,thisZ,thisT)=ratioStack(:,:,ratioAll(1)*thisZ,ratioAll(2)*thisT);
            end
        end
        stack=imresize(bigStack,1/ratioAll(1),'nearest');
    end
elseif dim==3
    if ratio(1)==1
%         stack=bwdist(stack);
        [x,y,z,t]=size(stack);
        for thisT=1:t
            stack(:,:,:,thisT)=bwdist(stack(:,:,:,thisT));
        end
    else
        if ndims(stack)==4
            [x,y,z,t]=size(stack);
            resStack=zeros(x,y,z,t);
            for thisT=1:t
                if ratio(1)>1
                    ratioStack=zeros(x,y,ratio(1)*z);
                    for thisZ=1:z
                        ratioStack(:,:,ratio(1)*thisZ)=stack(:,:,thisZ,thisT);
                    end
                    ratioStack=bwdist(ratioStack);
                    for thisZ=1:z
                        resStack(:,:,thisZ,thisT)=ratioStack(:,:,ratio(1)*thisZ);
                    end
                elseif ratio(1)<1
                    ratioStack=stack(:,:,:,thisT);
                    ratioStack=imresize(ratioStack,1/ratio(1),'nearest');
                    ratioStack=bwdist(ratioStack);
                    ratioStack=imresize(ratioStack,ratio(1),'nearest');
                    resStack(:,:,:,thisT)=ratioStack;
                end
            end
            stack=resStack;
        elseif ndims(stack)==3
            [x,y,z]=size(stack);
            if ratio(1)>1
                ratioStack=zeros(x,y,ratio(1)*z);
                for thisZ=1:z
                    ratioStack(:,:,ratio(1)*thisZ)=stack(:,:,thisZ);
                end
                ratioStack=bwdist(ratioStack);
                stack=zeros(x,y,z);
                for thisZ=1:z
                    stack(:,:,thisZ)=ratioStack(:,:,ratio(1)*thisZ);
                end
            elseif ratio(1)<1
                stack=imresize(stack,1/ratio(1),'nearest');
                stack=bwdist(stack);
                stack=imresize(stack,ratio(1),'nearest');

            end
        end
    end
end



end

