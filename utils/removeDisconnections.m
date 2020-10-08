%% actually removes all unconnected parts
function [image,rmFlag,connFlag] = removeDisconnections2(image)
%     disp('Removing disconnected flow partitions')
    rmFlag=false;
    connFlag=false;
    Nx=size(image,1);
    Ny=size(image,2);
    Nz=size(image,3);
%     geoSmaller(geoSmaller>1)=1;
    image=~image;
    CC=bwconncomp(image,6);
%     disp(['Identification complete: ',num2str(CC.NumObjects),' partitions identified'])
    % check each partition for end to end opening. if a partition is a
    %single part, and open on both sides, it must be connecte
    %throughout
%     disp('Filtering...')
    %% remove small blobs
    for i=1:CC.NumObjects
        if numel(CC.PixelIdxList{i})<Nz
%             disp(['Identified partition ',num2str(i),' as too short'])
            CC.PixelIdxList{i}=[];
            rmFlag=true;
%                 geoTemp=zeros(size(geoSmaller));
%                 geoTemp=geoTemp(:);
        end
    end  
    CC.PixelIdxList=CC.PixelIdxList(~cellfun('isempty',CC.PixelIdxList));
    %% check blob inlet-outlet hydraulic connectivity
    if numel(CC.PixelIdxList(:))>0
        for i=1:numel(CC.PixelIdxList(:))
            [~,~,z]=ind2sub([Nx Ny Nz],CC.PixelIdxList{i});
            numIn=sum(z==1);
            numOut=sum(z==Nz);
            if numIn==0
                CC.PixelIdxList{i}=[];
                rmFlag=true;
            elseif numOut==0
                CC.PixelIdxList{i}=[];
                rmFlag=true;
            end
            if rmFlag
%                 disp(['Identified large partition ',num2str(i),' as disconnected']) 
            end
        end
    end
    CC.PixelIdxList=CC.PixelIdxList(~cellfun('isempty',CC.PixelIdxList));
%     disp(['Removal complete: ',num2str(numel(CC.PixelIdxList(:))),' connected partitions identified'])
%     assert(numel(CC.PixelIdxList(:))>0,'Error: sample is disconnected')
    if numel(CC.PixelIdxList(:))>0
        connFlag=true;
    end
    CC.PixelIdxList=cell2mat(CC.PixelIdxList');
    image=~image;

    %%
    image=image(:);
    image(:)=true;
    image(CC.PixelIdxList)=false;
    
    image=reshape(image,[Nx Ny Nz]);
end