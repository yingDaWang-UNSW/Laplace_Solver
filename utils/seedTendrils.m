function [anchorPoints] = seedTendrils(temp)
            temp=double(~temp); %pores are 1
            midX=[];
            midY=[];
            CC=bwconncomp(temp,6);
%             disp(['Identification complete: ',num2str(CC.NumObjects),' partitions identified'])
            %check connectivity
            for j=1:numel(CC.PixelIdxList(:))
                rmFlag=false;
                [x,y,z]=ind2sub(size(temp),CC.PixelIdxList{j});
                numIn=sum(z==1);
                numOut=sum(z==size(temp,3));
                if numIn==0
                    rmFlag=true;
                elseif numOut==0
                    rmFlag=true;
                end
                if rmFlag %if disconnected, add tendrils at the median XY
                    midX=[midX;
                          round(quantile(x,0.5))];
                    midY=[midY;
                          round(quantile(y,0.5))];
                end
            end
            anchorPoints=[midX midY];
%             temp(midX,midY,:)=1;
%             [~,rmFlag] = removeDisconnections(temp==0);
%             assert(rmFlag==0,'Seeding failed, domain still disconnected')
end