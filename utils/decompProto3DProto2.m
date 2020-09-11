%% 3d version of domain decomp - unvectorised parameter extension - no onverlapping dual grid
%% things to improve:
% velocity fields should be calculated in the loop for partition flow rate,
% not recovered outside the loop.
% 

function [data] = decompProto3DProto2...
    (geopack,truncx,truncy,truncz,alpha,voxelSize,YDFlag,dFlag,dFile,dSaveFlag,...
    numPartitions,dualLength,numIterations,relaxFlag,relaxIndex,...
    precondFlag,pressureBasisFlag,upscalingFactor,microFlag,agmgFlag,...
    unixFlag,condDumpSerialiserFlag,diskDumpSerialiserFlag,preProcDiskDumpFlag,...
    diskFlag,runYDInSerial,runPartsInSerial,plotFlag)
    %% extra data for 3D
    data=struct();
    vtkFlag=false;
    resTol=0.0001;
    Nx=truncx(end);
    Ny=truncy(end);
    Nz=truncz(end);
    numPartitionsX=numPartitions(1);
    numPartitionsY=numPartitions(2);
    numPartitionsZ=numPartitions(3);
    incX=double(numPartitionsX>1);
    incY=double(numPartitionsY>1);
    incZ=double(numPartitionsZ>1);
    singlePassFlag=numPartitionsX+numPartitionsY+numPartitionsZ==3;
    dualLengthZ=dualLength;
    dualLengthY=ceil(max(truncy(end)-truncy(1)+1)/numPartitionsY/2); % allows for odd numbers. facility to safely overlap is not implemented
    dualLengthX=ceil(max(truncx(end)-truncx(1)+1)/numPartitionsX/2);
    %% prepare interim data file
    fileID = fopen('decompTempData.txt','a');
    fprintf(fileID,['\r\n\r\n Decompositional Algorithm 3D running for new test on ',geopack,...
                    ' sized ',num2str(max(truncx)),' ',num2str(max(truncy)),' ',num2str(max(truncz)) ,...
                    ' split into ',num2str(numPartitions),' parts with dual visibility of ',num2str(dualLength),' on ',datestr(datetime),'\r\n\r\n']);
    %% load binary
    tic
%     profile off
%     profile -memory on;
    %%
    disp('Loading Data')
    if dFlag
        geoSmaller=load([pwd,'/',geopack,'.mat'],'-mat');
        geoSmaller=geoSmaller.(geopack);

%         geoSmaller=load(['F:\Sandstone_ROI\',geopack,'.mat'],'ROIsandporesser');
%         geoSmaller=geoSmaller.ROIsandporesser;
%         geoSmaller=hlp_deserialize(geoSmaller);
%             geoSmaller=serialisedloadsave(geopack,[],false,false);
%             geoSmaller=geoSmaller.data;
%         geoSmaller=permute(geoSmaller,[2 1 3]);
        if truncx(end)+truncy(end)+truncz(end)-(size(geoSmaller))*ones(3,1)~=0
            geoSmaller=geoSmaller(truncx,truncy,truncz);
        end
    %     geoSmaller=permute(geoSmaller,[1 3 2]);
    % %%
%         for i=1:10:size(geoSmaller,1)
%             slice=squeeze(geoSmaller(:,:,i));
%             figure(1)
%             imagesc(slice);
%             drawnow
%         end
%         rock.perm=geoSmaller(:);
%         indexMap=[1:Nx*Ny*Nz]';
    end
    %% set secondary flags
    % clear geopack
    %voxelSize=10e-6;
    %dump during preconditioning to improve ram efficiency
%     precondDiskDumpFlag=true;
    if preProcDiskDumpFlag==true
        disp('Preprocessing Disk Dumping is Active')
    end
    if diskFlag
        preProcDiskDumpFlag=true;
    end
    anchorFlag=false;
    if agmgFlag
        if ~microFlag
            anchorFlag=true;
        end
    end
%     anchorFlag=false;
    if microFlag %force off disk dumping if using micropores [not necessary]
        preProcDiskDumpFlag=false;
        diskFlag=false;
    end
%     anchorThickness=1;
%     padThickness=1;
%     aggregateThickness=1;.
%     voxelsPerCoarseCell=250^3;
%     numPartitions = 5;%max(2,ceil(sum(geoSmaller(:))/voxelsPerCoarseCell));
%     maxZ=max(truncz)-mod(max(truncz),numPartitions);
%     truncz=truncz(1:maxZ);
%     geoSmaller=geoSmaller(truncx,truncy,truncz);
    %dualLength = ceil(max(truncz)/numPartitions/2);
    %microFlag=false;
    %agmgFlag=true;
    %precondFlag=false;
    %pressureBasisFlag=false;
%     RF=2;
    %relaxFlag=false;
%     maxRelax=max(numPartitionsZ/4,10);
%     minRelax=2;
%     upscaleRatio=1;
    %unixFlag=false;
    %runInSerial = true;
%     plotFlag=true;
%     plotFlagFig=false;
%     plotVelFlag=false;
    %diskFlag=true;
    if unixFlag
        rootDir='/media/user/HDD1/';
    else 
        rootDir='F:/';
    end
    diskDir=[rootDir,'transDump/'];
    % geoSmaller=zeros(10,10,10);
    % geoSmaller([1 10], :,:)=1;
    % geoSmaller(:,[1 10],:)=1;
    if microFlag
        delta=eps;
    else
        delta=0;
    end
    if runYDInSerial
      parDforArg = 0;
    else
      parDforArg = Inf;
    end
    if runPartsInSerial
      parForArg = 0;
    else
      parForArg = Inf;
    end
    numDualPartitionsZ=numPartitionsZ+incZ;
    numDualPartitionsY=numPartitionsY+incY;
    numDualPartitionsX=numPartitionsX+incX;
    %% remove cells
    if ~microFlag
        if dFlag
        [geoSmaller,~] = removeDisconnections(geoSmaller);
        end
    end
    clear CC indexMap temp slice z
    %% calculate weights
    if dFlag
        disp('Calculating fdgpa conductivity')
        D=bwdist(geoSmaller);
    %     alpha=1;
    %     delete(gcp('nocreate'))
    %     parpool('local',32)
        if YDFlag
            Dmax=calcDmax(D,parDforArg); 
        else
            maxdist=max(D(:));
            [Dmax] = FindWW2(D, maxdist,Nx, Ny, Nz, alpha,voxelSize );
        end
        D=alpha.*voxelSize.^2./8.*(2.*Dmax.*D-D.^2);
    else
        disp('Loading fdgpa conductivity')
%         D=load([dFile,'.mat']);
%         D=D.data;
        D=serialisedloadsave([rootDir,'condDump/',dFile,'.mat'],[],0,condDumpSerialiserFlag);
    end
    if dSaveFlag
        disp('Saving D to disk')
%         save([geopack,'3d.mat'],'D','-v7.3')
        serialisedloadsave([rootDir,'condDump/',geopack,'3d.mat'],D,1,condDumpSerialiserFlag);
    end
    if truncx(end)+truncy(end)+truncz(end)-(size(D))*ones(3,1)~=0
        D=D(truncx,truncy,truncz);
    end
        %     D=D(1:Nx,1:Ny,1:Nz);
%     Dmax=max(D(:)).*ones(size(D));
%     delete(gcp('nocreate'))
%     parpool('local',32)
%         for i=1:1:size(temp,3)
%             slice=squeeze(temp(:,:,i));
%             figure(1)
%             imagesc(slice);
%             drawnow
%         end
    clear Dmax
    %% generate intial guess from MS method
    if pressureBasisFlag
        disp('Basis estimation flag active, estimating pressure field from upscaled domain')
%         upscalingFactor=2;
        xDown=Nx/upscalingFactor;
        yDown=Ny/upscalingFactor;
        zDown=Nz/upscalingFactor;
        DDown=zeros(xDown,yDown,zDown);
        for i=1:xDown
            for j=1:yDown
                for k=1:zDown
                    x=upscalingFactor*i-(upscalingFactor-1):upscalingFactor*i;
                    y=upscalingFactor*j-(upscalingFactor-1):upscalingFactor*j;
                    z=upscalingFactor*k-(upscalingFactor-1):upscalingFactor*k;
                    locD=D(x,y,z);
        %             locD(:)
                    %% arithmetic average conductivity
                    condArithLoc=mean(locD(:));
                    %% Agarwal18 method
        %             condAgarLoc=
                    %% other method?

                    %% distribute to DDown
                    DDown(i,j,k)=condArithLoc;
                end
            end
        end
        % generate grid
        geoDownCustom=DDown==0;
        disp('Seeding domain with anchors')  
        anchors=false(size(DDown));
        disp(['Seeding Tendrils'])
        [anchorPoints] = seedTendrils(geoDownCustom);
        if ~isempty(anchorPoints)
            for m=1:size(anchorPoints,1)
                x=anchorPoints(m,1);
                y=anchorPoints(m,2);
        %                                     anchors(x,y,cellCellIndsZ{k})=true;
                anchors(x,y,:)=true;
            end
        end
        disp('Superimposing anchors with weights')
        DDown(anchors)=DDown(anchors)+eps;
        geoDownCustom=DDown==0;
        [KDownCustom,solDown,GDown]= cond4trai(geoDownCustom,voxelSize*upscalingFactor,1,false,true,false,false,0,1,true,DDown);
        upscaleFaceCentroids=GDown.faces.centroids;
        upscaleFacePressures=solDown.facePressure;
        upscaleFaceIndexMap=(1:GDown.faces.num)';
        clear GDown solDown
        
    end
    
    %% anchor image with eps pseudo connectivity
    anchors=false(size(D));
    if ~plotFlag
        geoSmaller=D==0;
    end
    intervalX=cell(numPartitionsX,1);
    intervalY=cell(numPartitionsY,1);
    intervalZ=cell(numPartitionsZ,1);
    for i=1:numPartitionsX
        intervalX{i}=i*Nx/numPartitionsX-Nx/numPartitionsX+1:i*Nx/numPartitionsX;
    end
    for j=1:numPartitionsY
        intervalY{j}=j*Ny/numPartitionsY-Ny/numPartitionsY+1:j*Ny/numPartitionsY;
    end
    for k=1:numPartitionsZ
        intervalZ{k}=k*Nz/numPartitionsZ-Nz/numPartitionsZ+1:k*Nz/numPartitionsZ;
    end

    dualIntervalX=cell(numDualPartitionsX,1);
    dualIntervalY=cell(numDualPartitionsZ,1);
    dualIntervalZ=cell(numDualPartitionsZ,1);
    iterNumX=1;
    iterNumY=1;
    iterNumZ=1;
    visibilityIntervalX=dualLengthX;
    visibilityIntervalY=dualLengthY;
    visibilityIntervalZ=dualLengthZ;

    for i=1:numDualPartitionsX
        if i==1
            if numDualPartitionsX==1
                dualIntervalX{i}=iterNumX:Nx;
            else
                dualIntervalX{i}=iterNumX:visibilityIntervalX;
            end
        elseif i==numDualPartitionsX
            dualIntervalX{i}=max(intervalX{i-1})-visibilityIntervalX+1:max(intervalX{i-1});
        else
            dualIntervalX{i}=min(intervalX{i})-visibilityIntervalX:min(intervalX{i})+visibilityIntervalX-1;
        end
        iterNumX=max(dualIntervalX{i})+1;
    end

    for i=1:numDualPartitionsY
        if i==1
            if numDualPartitionsY==1
                dualIntervalY{i}=iterNumY:Ny;
            else
                dualIntervalY{i}=iterNumY:visibilityIntervalY;
            end
        elseif i==numDualPartitionsY
            dualIntervalY{i}=max(intervalY{i-1})-visibilityIntervalY+1:max(intervalY{i-1});
        else
            dualIntervalY{i}=min(intervalY{i})-visibilityIntervalY:min(intervalY{i})+visibilityIntervalY-1;
        end
        iterNumY=max(dualIntervalY{i})+1;
    end

    for i=1:numDualPartitionsZ
        if i==1
            if numDualPartitionsZ==1
                dualIntervalZ{i}=iterNumZ:Nz;
            else
                dualIntervalZ{i}=iterNumZ:visibilityIntervalZ;
            end
        elseif i==numDualPartitionsZ
            dualIntervalZ{i}=max(intervalZ{i-1})-visibilityIntervalZ+1:max(intervalZ{i-1});
        else
            dualIntervalZ{i}=min(intervalZ{i})-visibilityIntervalZ:min(intervalZ{i})+visibilityIntervalZ-1;
        end
        iterNumZ=max(dualIntervalZ{i})+1;
    end
                
                
    if ~microFlag
        if agmgFlag
            if anchorFlag
                disp('Seeding domain with anchors')       
                for i=1:numPartitionsX
                    for j=1:numPartitionsY
                        for k=1:numPartitionsZ
                            temp=geoSmaller(intervalX{i},intervalY{j},intervalZ{k});
                            localOrigin=[min(intervalX{i})-1,min(intervalY{j})-1,min(intervalZ{k})-1];
                            disp(['Seeding Tendrils for grid block',num2str(i),num2str(j),num2str(k)])
                            [anchorPoints] = seedTendrils(temp);
                            if ~isempty(anchorPoints)
                                for m=1:size(anchorPoints,1)
                                    x=anchorPoints(m,1)+localOrigin(1);
                                    y=anchorPoints(m,2)+localOrigin(2);
%                                     anchors(x,y,cellCellIndsZ{k})=true;
                                    anchors(x,y,:)=true;
                                end
                            end
                        end
                    end
                end
        %         rock.perm=anchors(:);
        %         indexMap=[1:Nx*Ny*Nz]';
        %     %     %needed due to obnoxious memory redundancy associated with G
        %         G=generateReducedGridYDW(size(geoSmaller),indexMap(rock.perm==0|isnan(rock.perm)));
        %         plotGrid(G)
               for i=1:numDualPartitionsX
                    for j=1:numDualPartitionsY
                        for k=1:numDualPartitionsZ
                            temp=geoSmaller(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
                            localOrigin=[min(dualIntervalX{i})-1,min(dualIntervalY{j})-1,min(dualIntervalZ{k})-1];
                            disp(['Seeding Tendrils for dual grid block',num2str(i),num2str(j),num2str(k)])
                            [anchorPoints] = seedTendrils(temp);
                            if ~isempty(anchorPoints)
                                for m=1:size(anchorPoints,1)
                                    x=anchorPoints(m,1)+localOrigin(1);
                                    y=anchorPoints(m,2)+localOrigin(2);
%                                     anchors(x,y,intervalZ{k})=true;
                                    anchors(x,y,:)=true;
                                end
                            end
                        end
                    end
                end       
            end
%             rock.perm=anchors(:);
%             indexMap=[1:Nx*Ny*Nz]';
%             %needed due to obnoxious memory redundancy associated with G
%             G=generateReducedGridYDW(size(geoSmaller),indexMap(rock.perm==0|isnan(rock.perm)));
%             figure(60);
%             plotGrid(G)
            disp('Superimposing anchors with weights')
            D(anchors)=D(anchors)+eps;
            if plotFlag
                dAnchors=D.*-anchors;
                dAnchors=dAnchors<0;
            end
        end
    end
    clear temp

    %% split weights and binary and NDCF arrays
    disp('Splitting coarse partitions')
    % geoSmallerParts=cell(numPartitions,1);
%     Dparts=cell(numPartitions,1);
    Gparts=cell(numPartitionsX,numPartitionsY,numPartitionsZ);
    transParts=cell(numPartitionsX,numPartitionsY,numPartitionsZ);
    faceLocsX=zeros(numPartitionsX,2);
    faceLocsY=zeros(numPartitionsY,2);
    faceLocsZ=zeros(numPartitionsZ,2);
    dirPressVecX=linspace(1,0,Nz)';
    dirPressVecY=linspace(1,0,Nz)';
    dirPressVecZ=1;
    %parallel
    for i=1:numPartitionsX
        faceLocsX(i,:)=[min(intervalX{i})-1 max(intervalX{i})];
%         dirPressVecX=[dirPressVecX linspace(1,0,Nx)'];
    end
    for i=1:numPartitionsY
        faceLocsY(i,:)=[min(intervalY{i})-1 max(intervalY{i})];
%         dirPressVecY=[dirPressVecY linspace(1,0,Ny)'];
    end
    % perpendicular
    for i=1:numPartitionsZ
        faceLocsZ(i,:)=[min(intervalZ{i})-1 max(intervalZ{i})];
        dirPressVecZ=[dirPressVecZ (1-0)*(numPartitionsZ-i)/numPartitionsZ];
    end
    
    
    %% map centroids from upscaled to fine scale
    if pressureBasisFlag && ~singlePassFlag
        disp("Assigning Pressure Basis to Initial Domain")
        %% get the upscaled faces of interest 
        upscaleFaceCentroids=round(upscaleFaceCentroids,1);
        FOIX=unique(faceLocsX(:))./upscalingFactor;
        FOIX=FOIX(2:end-1);
        FOIY=unique(faceLocsY(:))./upscalingFactor;
        FOIY=FOIY(2:end-1);
        FOIZ=unique(faceLocsZ(:))./upscalingFactor;
        FOIZ=FOIZ(2:end-1);
        internalFaceCentroids=[];
        internalFacePressures=[];
        for i=1:numel(FOIX)
            inds=upscaleFaceCentroids(:,1)==FOIX(i);
            internalFaceCentroids=[internalFaceCentroids; upscaleFaceCentroids(inds,:)];
            internalFacePressures=[internalFacePressures; upscaleFacePressures(inds)];
%             figure(1);plotFaces(GDown,upscaleFaceIndexMap(inds));hold on
        end
        for i=1:numel(FOIY)
            inds=upscaleFaceCentroids(:,2)==FOIY(i);
            internalFaceCentroids=[internalFaceCentroids; upscaleFaceCentroids(inds,:)];
            internalFacePressures=[internalFacePressures; upscaleFacePressures(inds)];

%             figure(1);plotFaces(GDown,upscaleFaceIndexMap(inds));hold on
        end
        for i=1:numel(FOIZ)
            inds=upscaleFaceCentroids(:,3)==FOIZ(i);
            internalFaceCentroids=[internalFaceCentroids; upscaleFaceCentroids(inds,:)];
            internalFacePressures=[internalFacePressures; upscaleFacePressures(inds)];

%             figure(1);plotFaces(GDown,upscaleFaceIndexMap(inds));hold on
        end
    end
    %%
    gridNDFP=cell(numPartitionsX,numPartitionsY,numPartitionsZ);
    rmFlags=zeros(numPartitionsX,numPartitionsY,numPartitionsZ);
    for i=1:numPartitionsX
        for j=1:numPartitionsY
            for k=1:numPartitionsZ
                rock=struct();
                disp(['Splitting coarse partition ', num2str(i),num2str(j),num2str(k)])
                disp('Generating representative coarse grid')
                temp=D(intervalX{i},intervalY{j},intervalZ{k});
                if microFlag
                    if i+j+k==3
                        disp('Microporosity flag active, extracting single grid')
                        G=cartGrid(size(double(temp)));
                        G=mcomputeGeometry(G);
                    end
                    transParts{i,j,k}=computeTrans(G,struct('perm',double(temp(:))+delta),'verbose',true); 
                else
                    G=[];
                    disp('Microporosity inactive, generating reduced cartesian grid set for each partition')
                    [~,rmFlag] = removeDisconnections(~temp);
                    if rmFlag==1
                            pause(1)
                    end
                    rmFlags(i,j,k)=rmFlag;
                    rock.perm=double(temp(:))+delta;
                    indexMap=(1:numel(rock.perm))';
                    Gparts{i,j,k}=generateReducedGridYDW(size(temp),indexMap(rock.perm==0|isnan(rock.perm)));
                    rock.perm=rock.perm(rock.perm>0);
                    Gparts{i,j,k}=mcomputeGeometry(Gparts{i,j,k});
                    if plotFlag
                        if agmgFlag
                        Gparts{i,j,k}.anchors=dAnchors(intervalX{i},intervalY{j},intervalZ{k});
                        Gparts{i,j,k}.anchors=Gparts{i,j,k}.anchors(:);
                        Gparts{i,j,k}.anchors=Gparts{i,j,k}.anchors(Gparts{i,j,k}.cells.indexMap);
                        end
                    end
                    transParts{i,j,k}=computeTrans(Gparts{i,j,k},rock,'verbose',true);
                    %% find pressure mapping 
                    if pressureBasisFlag && ~singlePassFlag
                        disp(['Mapping Basis Pressures to block ',num2str(i),num2str(j),num2str(k)])
                        localOrigin=[faceLocsX(i,1),faceLocsY(j,1),faceLocsZ(k,1)];
                        primalFaceCentroids=Gparts{i,j,k}.faces.centroids+localOrigin; 
                        % for each of the 6 faces, find the upscaled face
                        % values in the order [xmin xmax; ymin ymax; zmin zmax]
                        boundingBox=[faceLocsX(i,:);
                                     faceLocsY(j,:);
                                     faceLocsZ(k,:)];

                        inactiveFaces=[i==1,i==numPartitionsX;
                                       j==1,j==numPartitionsY;
                                       k==1,k==numPartitionsZ];
                        NDFP=cell(3,2);
                        tic
                        for m=1:3 %find primal faces, then map each face to the upscaled value [no left and right nonsense]
                            for n=1:2
                                otherDims=setdiff([1,2,3],m);
    %                             otherBound=setdiff([1,2],n);
                                if inactiveFaces(m,n)==1
                                    NDFP{m,n}=[];
                                    continue
                                end
                                primalFaceInds=primalFaceCentroids(:,m)==boundingBox(m,n); %retains ordering built in G

    %                             subFaceCentroids=G.faces.centroids(faceIndexList(G.faces.centroids(:,m)==bc(m,n)),:);
                                tempFaceCentroids=primalFaceCentroids(primalFaceInds,:)./upscalingFactor; 
                                dirMask=mod(tempFaceCentroids,1)~=0;
                                tempFaceCentroids(dirMask)=floor(tempFaceCentroids(dirMask))+0.5;
                                tempFaceCentroids=round(tempFaceCentroids,1);
                                tempPressureVec=zeros(size(dirMask,1),1);
                                for x=1:size(dirMask,1)
%                                     disp(['Progress: ',num2str(x/numel(dirMask))])
    %                                 [~,lookupInd]=ismember(tempFaceCentroids(x,:), upscaleFaceCentroids,'rows'); %this is so fucking slow. use dimensional squeeze.
                                    lookupInd=upscaleFaceCentroids(:,1)==tempFaceCentroids(x,1) &...% squeeze along x dim
                                              upscaleFaceCentroids(:,2)==tempFaceCentroids(x,2) &...; % bound along other dims
                                              upscaleFaceCentroids(:,3)==tempFaceCentroids(x,3); %works for doubles moraculously!
                                    dirPressVal=tempFaceCentroids(x,3)/zDown;
                                    dirPressWeight=0;
                                    try%if ~isempty(upscaleFacePressures(lookupInd))
                                        tempPressureVec(x,1)=((1-dirPressWeight)*upscaleFacePressures(lookupInd)+dirPressWeight*dirPressVal);
                                    catch
                                        tempPressureVec(x,1)=0.5;
                                    end
                                end
                                NDFP{m,n}=[primalFaceCentroids(primalFaceInds,:), tempPressureVec];
                                %from the upscaled centroids, find the faces of
                                %interest and sort the values in the same order
                                %as the primal fcae indexes
    %                             temp=zeros(numel(faceIndexMap),1);
    %                             temp(primalFaceInds)=tempPressureVec;
%                                 tempG=Gparts{i,j,k};
%                                 tempG.nodes.coords=tempG.nodes.coords+localOrigin;
    %                             figure(1);plotFaces(tempG,faceIndexMap(primalFaceInds),'FaceVertexCData',tempPressureVec,'faceColor','flat');hold on;colorbar;caxis([0,1]); axis tight; axis equal;
    %                             upscaleFaceInds
                            end
                        end
                        toc
                        gridNDFP{i,j,k}=NDFP;
                    end
                    if preProcDiskDumpFlag
                        disp('Disk dumping active, dumping trans and grid data to disk')
            %             temp=transParts{i};
            %             save([diskDir,'trans',num2str(i),'.mat'],'temp','-v7.3')
                        serialisedloadsave([diskDir,'trans', num2str(i),num2str(j),num2str(k),'.mat'],transParts{i,j,k},1,diskDumpSerialiserFlag);
                        transParts{i,j,k}=[];
                        serialisedloadsave([diskDir,'grid', num2str(i),num2str(j),num2str(k),'.mat'],Gparts{i,j,k},1,diskDumpSerialiserFlag);
                        Gparts{i,j,k}=[];
                    end
                end
            end
        end
    end
    if pressureBasisFlag && ~singlePassFlag
      
    else
        downP=[];
        gridNDFP=[];
    end
    clear Dparts temp temp2 rock 
    %% split into dual arrays
    disp('Splitting dual coarse partitions')
    Gdualparts=cell(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ);
    transDualParts=cell(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ);
    dualFaceLocsX=zeros(numDualPartitionsX,2);
    dualFaceLocsY=zeros(numDualPartitionsY,2);
    dualFaceLocsZ=zeros(numDualPartitionsZ,2);
    dirDualPressVecZ=zeros(numDualPartitionsZ,2);
    %parallel
    for i=1:numDualPartitionsX
        dualFaceLocsX(i,:)=[min(dualIntervalX{i})-1 max(dualIntervalX{i})];
    end
    for i=1:numDualPartitionsY
        dualFaceLocsY(i,:)=[min(dualIntervalY{i})-1 max(dualIntervalY{i})];
    end
    % perpendicular
    for i=1:numDualPartitionsZ
        dualFaceLocsZ(i,:)=[min(dualIntervalZ{i})-1 max(dualIntervalZ{i})];
        dirDualPressVecZ(i,:)=[(1-0)*(Nz-min(dualIntervalZ{i})+1)/Nz (1-0)*(Nz-max(dualIntervalZ{i}))/Nz];
    end
    rmDualFlags=zeros(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ);
    if ~singlePassFlag
        for i=1:numDualPartitionsX
            for j=1:numDualPartitionsY
                for k=1:numDualPartitionsZ
                    rock=struct();
                    disp(['Splitting coarse dual partition ', num2str(i),num2str(j),num2str(k)])
                    disp('Generating representative dual coarse grid')
                    temp=D(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
                    if microFlag
                        disp('Microporosity flag active, generating only trans and repeated grids')
                        %find corner, edge, face, and body grids but also
                        %remove based on partitioning dimensions
                        dimVec=[i==1||i==numDualPartitionsX,...
                                j==1||j==numDualPartitionsY,...
                                k==1||k==numDualPartitionsZ];
                        dimMask=[numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ];
                        dimMask=dimMask>1;
                        dimVec=dimVec.*dimMask;
                        if sum(dimMask)==2
                            G8=[];
                        elseif sum(dimMask)==1
                            G4=[];
                            G8=[];
                        end
                        % handle eighths [3 dims must be 1 or end] 
                        if sum(dimVec)==3
                            if ~exist('G8','var')
                                disp('Microporosity flag active, extracting single eighth grid')
                                G8=cartGrid(size(double(temp)));
                                G8=mcomputeGeometry(G8);
                            end
                            transDualParts{i,j,k}=computeTrans(G8,struct('perm',double(temp(:))+delta),'verbose',true); 
                        %handle fulls [no dims must be 1 or end]
                        elseif sum(dimVec)==0
                            if ~exist('G1','var')
                                disp('Microporosity flag active, extracting single full grid')
                                G1=cartGrid(size(double(temp)));
                                G1=mcomputeGeometry(G1);
                            end
                            transDualParts{i,j,k}=computeTrans(G1,struct('perm',double(temp(:))+delta),'verbose',true); 
                        %handle halves [only 1 dim must be 1 or end]
                        elseif sum(dimVec)==1
                            gridSize=size(temp);
                            [~,gridInd]=max(dimVec);
                            if ~exist('GH','var')
                                GH=cell(3,1);
                                disp('Microporosity flag active, extracting triple half grid')
                                GH{1}=cartGrid(gridSize.*(dimVec+1)./[2,1,1]);
                                GH{1}=mcomputeGeometry(GH{1});
                                GH{2}=cartGrid(gridSize.*(dimVec+1)./[1,2,1]);
                                GH{2}=mcomputeGeometry(GH{2});
                                GH{3}=cartGrid(gridSize.*(dimVec+1)./[1,1,2]);
                                GH{3}=mcomputeGeometry(GH{3});
                            end
                            transDualParts{i,j,k}=computeTrans(GH{gridInd},struct('perm',double(temp(:))+delta),'verbose',true); 
                        %handle quarters [else]
                        elseif sum(dimVec)==2
                            gridSize=size(temp);
                            [~,gridInd]=min(dimVec);
                            if ~exist('G4','var')
                                G4=cell(3,1);
                                disp('Microporosity flag active, extracting triple quarter grid')
%                                 G4=cartGrid(size(double(temp)));
%                                 G4=mcomputeGeometry(G4); 
                                G4{1}=cartGrid(gridSize.*(dimVec+1)./[1,2,2]);
                                G4{1}=mcomputeGeometry(G4{1});
                                G4{2}=cartGrid(gridSize.*(dimVec+1)./[2,1,2]);
                                G4{2}=mcomputeGeometry(G4{2});
                                G4{3}=cartGrid(gridSize.*(dimVec+1)./[2,2,1]);
                                G4{3}=mcomputeGeometry(G4{3});
                            end
                            transDualParts{i,j,k}=computeTrans(G4{gridInd},struct('perm',double(temp(:))+delta),'verbose',true); 
                        end
                    else
                        G1=[];GH=[];G4=[];G8=[];
                        disp('Microporosity inactive, generating reduced cartesian grid set for each dual partition')
                        [~,rmDualFlag] = removeDisconnections(~temp);
                        if rmDualFlag==1
                                pause(1)
                        end
                        rmDualFlags(i,j,k)=rmDualFlag;
                        rock.perm=double(temp(:))+delta;
                        indexMap=(1:numel(rock.perm))';
                        Gdualparts{i,j,k}=generateReducedGridYDW(size(temp),indexMap(rock.perm==0|isnan(rock.perm)));
                        rock.perm=rock.perm(rock.perm>0);
                        Gdualparts{i,j,k}=mcomputeGeometry(Gdualparts{i,j,k});
                        if plotFlag
                            if agmgFlag
                            Gdualparts{i,j,k}.anchors=dAnchors(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
                            Gdualparts{i,j,k}.anchors=Gdualparts{i,j,k}.anchors(:);
                            Gdualparts{i,j,k}.anchors=Gdualparts{i,j,k}.anchors(Gdualparts{i,j,k}.cells.indexMap);
                            end
                        end
                        transDualParts{i,j,k}=computeTrans(Gdualparts{i,j,k},rock,'verbose',true); 
                        if preProcDiskDumpFlag
                            disp('Disk dumping active, dumping trans and grid data to disk')
                %             temp=transParts{i};
                %             save([diskDir,'trans',num2str(i),'.mat'],'temp','-v7.3')
                            serialisedloadsave([diskDir,'dualTrans', num2str(i),num2str(j),num2str(k),'.mat'],transDualParts{i,j,k},1,diskDumpSerialiserFlag);
                            transDualParts{i,j,k}=[];
                            serialisedloadsave([diskDir,'dualGrid', num2str(i),num2str(j),num2str(k),'.mat'],Gdualparts{i,j,k},1,diskDumpSerialiserFlag);
                            Gdualparts{i,j,k}=[];
                        end
                    end
                end
            end
        end
    end
    %% face locations
    disp('Setting face locations') % regroup boundary face locations to suit grids
    if numDualPartitionsX>1
        cellDualFaceLocsX=dualFaceLocsX(2:end-1);
        cellDualFaceLocsX=reshape(cellDualFaceLocsX,[numPartitionsX 2]);
        cellDualFaceLocsX=fliplr(cellDualFaceLocsX);
    else
        cellDualFaceLocsX=dualFaceLocsX;
    end

    
    if numDualPartitionsY>1
        cellDualFaceLocsY=dualFaceLocsY(2:end-1);
        cellDualFaceLocsY=reshape(cellDualFaceLocsY,[numPartitionsY 2]);
        cellDualFaceLocsY=fliplr(cellDualFaceLocsY);
    else
        cellDualFaceLocsY=dualFaceLocsY;
    end

    
    if numDualPartitionsZ>1
        cellDualFaceLocsZ=dualFaceLocsZ(2:end-1);
        cellDualFaceLocsZ=reshape(cellDualFaceLocsZ,[numPartitionsZ 2]);
        cellDualFaceLocsZ=fliplr(cellDualFaceLocsZ);
    else
        cellDualFaceLocsZ=dualFaceLocsZ;
    end

    
    bcZ=faceLocsZ(1,:); % extract local boundary face locations
    bcX=faceLocsX(1,:); 
    bcY=faceLocsY(1,:); 
    
    bc=[bcX;bcY;bcZ];
    
%     dualhalfbcZ=dualFaceLocsZ(1,:); % extract local bc half dual faces
%     dualhalfbcX=dualFaceLocsX(1,:); 
%     dualhalfbcY=dualFaceLocsY(1,:); 
%     dualbcZ=dualhalfbcZ.*2; % extract local bc dual faces
%     dualbcX=dualhalfbcX.*2; % extract local bc dual faces
%     dualbcY=dualhalfbcY.*2; % extract local bc dual faces

    clear matPressure upscalePressure
%% reload grid and trans data if flags
    if preProcDiskDumpFlag==true
        if diskFlag==false
            disp('Loading all partition data back to RAM')
            for i=1:numPartitionsX
                for j=1:numPartitionsY
                    for k=1:numPartitionsZ
                        Gparts{i,j,k}=serialisedloadsave([diskDir,'grid', num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        transParts{i,j,k}=serialisedloadsave([diskDir,'trans', num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        if ~diskDumpSerialiserFlag
                            Gparts{i,j,k}=Gparts{i,j,k}.data;
                            transParts{i,j,k}=transParts{i,j,k}.data;
                        end
                        disp(['Partition ', num2str(i),num2str(j),num2str(k),' loaded back into RAM'])
                    end
                end
            end
            for i=1:numDualPartitionsX
                for j=1:numDualPartitionsY
                    for k=1:numDualPartitionsZ
                        Gdualparts{i,j,k}=serialisedloadsave([diskDir,'dualGrid',num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        transDualParts{i,j,k}=serialisedloadsave([diskDir,'dualTrans',num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        if ~diskDumpSerialiserFlag
                            Gdualparts{i,j,k}=Gdualparts{i,j,k}.data;
                            transDualParts{i,j,k}=transDualParts{i,j,k}.data;
                        end
                        disp(['Dual partition ',num2str(i),num2str(j),num2str(k),' loaded back into RAM'])
                    end
                end
            end
        end
    end
 
    %% exploded plotting routine
%     sepLen=0;
%     for i=1:numPartitionsX
%         for j=1:numPartitionsY
%             for k=1:numPartitionsZ
%                 originX=(i-1)*Nx/numPartitionsX+(i-1)*sepLen;
%                 originY=(j-1)*Nx/numPartitionsY+(j-1)*sepLen;
%                 originZ=(k-1)*Nx/numPartitionsZ+(k-1)*sepLen;
%                 tempG=Gparts{i,j,k};
%                 tempG.nodes.coords=tempG.nodes.coords+[originX originY originZ];
%                 figure(69)
%                 plotGrid(tempG);hold on;axis equal;axis tight
%             end
%         end
%     end
%     %% exploded plotting routine
%     sepLen=35;
%     for i=1:numDualPartitionsX
%         for j=1:numDualPartitionsY
%             for k=1:numDualPartitionsZ
%                 if i==1
%                     fixx=Nx/numDualPartitionsX/2+sepLen/8;
%                 else 
%                     fixx=0;
%                 end
%                 if j==1
%                     fixy=Ny/numDualPartitionsY/2+sepLen/8;
%                 else
%                     fixy=0;
%                 end
%                 if k==1
%                     fixz=Nz/numDualPartitionsZ/2+sepLen/8;
%                 else
%                     fixz=0;
%                 end
%                 originX=(i-1)*Nx/numDualPartitionsX+(i-1)*sepLen+fixx;
%                 originY=(j-1)*Nx/numDualPartitionsY+(j-1)*sepLen+fixy;
%                 originZ=(k-1)*Nx/numDualPartitionsZ+(k-1)*sepLen+fixz;
%                 tempG=Gdualparts{i,j,k};
%                 tempG.nodes.coords=tempG.nodes.coords+[originX originY originZ];
%                 figure(70)
%                 plotGrid(tempG);hold on;axis equal; axis tight
%             end
%         end
%     end
    %%  solver routine
    % iterate until matching!
    % generate coarse grids
    disp('Solving...')
    residual=100;
    iterNum=0;
    resVec=[];
    permVec=[];
    minPermVec=[];
    midqPermVec=[];
    balancePermVec=[];
    timeVec=[];
    qMeanVec=[];
    qBalVec=[];
    qMidVec=[];
    clear rock indexMap 
    fclose(fileID);
    simTime=tic;
    numTotParts=numPartitionsX*numPartitionsY*numPartitionsZ;
    numTotDualParts=numDualPartitionsX*numDualPartitionsY*numDualPartitionsZ;
    Area = (voxelSize)^2;

    while iterNum<=numIterations && residual>resTol
        fileID = fopen('decompTempData.txt','a');
        midBCPressures=cell(numTotParts,1);
        subVel=cell(numTotParts,1);
        subVelZ=cell(numTotParts,1);
        subVelX=cell(numTotParts,1);
        subVelY=cell(numTotParts,1);
        sumVelZ=cell(numTotParts,1);
        %solve flow with coarse grids
        iterNum=iterNum+1;
        disp('Solving fine partitions')
%         toc
%         tic
        localDual=[cellDualFaceLocsX(1,:);
                   cellDualFaceLocsY(1,:);
                   cellDualFaceLocsZ(1,:)];
        clear totVel
%         parfor (i=1:numPartitions,parforArg)%(idx = 1:N, parforArg)
%         for (i=1:numPartitions)%(idx = 1:N, parforArg)
        %% grid routine
%         parfor (n=1:numTotParts,parForArg)
        for (n=1:numTotParts)
            partTime=tic;
            [i,j,k]=ind2sub([numPartitionsX,numPartitionsY,numPartitionsZ],n);
%         for i=1:numPartitionsX
%             for j=1:numPartitionsY
%                 for k=1:numPartitionsZ
                    if rmFlags(i,j,k)==0
                        solver='AGMG';
                    else
                        solver='MLDIVIDE';
                    end
                    localOrigin=[faceLocsX(i,1),faceLocsY(j,1),faceLocsZ(k,1)];
                    dualFlag=false;
                    if ~diskFlag
                        S=transParts{i,j,k};
                        if ~microFlag
                            G2=Gparts{i,j,k};
                        else
                            G2=G;
                        end
                    else
                        disp(['Loading grid ',num2str(i),num2str(j),num2str(k),' to RAM'])
                        S=serialisedloadsave([diskDir,'trans',num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        G2=serialisedloadsave([diskDir,'grid',num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        if ~diskDumpSerialiserFlag
                            G2=G2.data;
                            S=S.data;
                        end
                    end
                    if iterNum==1 % populate NDFP with initial guess data
                        if ~pressureBasisFlag
                            NDFP{1,1}=dirPressVecX(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{1,2}=dirPressVecX(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{2,1}=dirPressVecY(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{2,2}=dirPressVecY(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{3,1}=dirPressVecZ(k);
                            NDFP{3,2}=dirPressVecZ(k+1); 
                        else
                            NDFP=gridNDFP{i,j,k};
                        end
                        init=[];
                    else
                        NDFP=gridNDFP{i,j,k};
                        init=[];
                    end
                    % determine active faces
                    activeFaces=[i==1,i==numPartitionsX;
                                 j==1,j==numPartitionsY;
                                 k==1,k==numPartitionsZ];
                    precond=(localOrigin(3)+Nz/numPartitionsZ/2)/2;
                    [state,dualPressure] =nonDirichletCorrectionSubroutine3D...
                        (G2,S,NDFP,bc,localOrigin,localDual,plotFlag,iterNum,...
                        dualFlag,rmFlags(i,j,k),~activeFaces,init,[]);
                    midBCPressures{n}=dualPressure;
                    disp(['Solved fine partition ',num2str(i),num2str(j),num2str(k),...
                          ' with ', solver,' ',num2str(toc(partTime)),': seconds elapsed'])
        %             disp('Extracting Flux')
                    tempVel=faceFlux2cellVelocity(G2,state.flux);
                    tempVel=tempVel(:,3);
                    cartVelZ=zeros(G2.cartDims);
                    cartVelZ=cartVelZ(:);
                    cartVelZ(G2.cells.indexMap)=tempVel;
                    tempVel=reshape(cartVelZ,G2.cartDims);
                    tempuz=zeros(1,1,G2.cartDims(3));
                    for (p = 1:G2.cartDims(3))
                        tempuz(1,1,p) = sum(sum(tempVel(:,:,p)));
                    end
                    sumVelZ{n}=tempuz;
                    if vtkFlag
                        subVel{n}=faceFlux2cellVelocity(G2,state.flux);
                        tempVel=subVel{n}(:,3);
                        cartVelZ=zeros(G2.cartDims);
                        cartVelZ=cartVelZ(:);
                        cartVelZ(G2.cells.indexMap)=tempVel;
                        subVelZ{n}=reshape(cartVelZ,G2.cartDims);
                        tempVel=subVel{n}(:,1);
                        cartVelX=zeros(G2.cartDims);
                        cartVelX=cartVelX(:);
                        cartVelX(G2.cells.indexMap)=tempVel;
                        subVelX{n}=reshape(cartVelX,G2.cartDims);

                        tempVel=subVel{n}(:,2);
                        cartVelY=zeros(G2.cartDims);
                        cartVelY=cartVelY(:);
                        cartVelY(G2.cells.indexMap)=tempVel;
                        subVelY{n}=reshape(cartVelY,G2.cartDims);
                    end
%                     if plotFlag
%                         totVel{i,j,k}=subVel{i,j,k};
%                     end

%                 end
%             end
%         end
        end
        clear G2
%         toc
%         tic
        disp('Reconstructing fluxes')
        %indexmapping will fix this here
        if vtkFlag
            subVelZ=reshape(subVelZ,[numPartitionsX,numPartitionsY,numPartitionsZ]);
            subVelZ=cell2mat(subVelZ);
            subVelZ(isnan(subVelZ)) = 0;
            subVelZ  = (1/voxelSize)*subVelZ;
            subVelY=reshape(subVelY,[numPartitionsX,numPartitionsY,numPartitionsZ]);
            subVelY=cell2mat(subVelY);
            subVelY(isnan(subVelY)) = 0;
            subVelY  = (1/voxelSize)*subVelY;

            subVelX=reshape(subVelX,[numPartitionsX,numPartitionsY,numPartitionsZ]);
            subVelX=cell2mat(subVelX);
            subVelX(isnan(subVelX)) = 0;
            subVelX  = (1/voxelSize)*subVelX;

            [X,Y,Z]=meshgrid(truncx,truncy,truncz);        
            vtkwrite(['tempVelsIter',sprintf('%04d',iterNum),'.vtk'],'structured_grid',X,Y,Z,'vectors','Velocity',subVelX,subVelY,subVelZ);
        end

%         subVel=reshape(subVel,[Nx Ny Nz]);
        tempuz = zeros(Nz,1);
%         parfor (k = 1:Nz,parforArg)
        sumVelZ=reshape(sumVelZ,[numPartitionsX,numPartitionsY,numPartitionsZ]);
        sumVelZ=cell2mat(sumVelZ);
        for (k = 1:Nz)
            tempuz(k) = sum(sum(sumVelZ(:,:,k)));
        end
        q = tempuz(:).*Area/voxelSize;
        q=q';
        meanq = mean(q);
        figure(1);clf(1);
        plot(1:numel(q),q);hold on;
        plot(1:numel(q),repmat(meanq,[numel(q),1]))
        data.q{iterNum,1}=q;
        plot(1:numel(q),data.q{1,1});
        %calculate mass imbalance
        fq=q-min(q);
        deltaq=max(q)-min(q);
        areaq=trapz(fq);
        for a=1:100
            tempq=fq-deltaq*a/100;
            tempqarea=trapz(tempq);
            if tempqarea<areaq/2
                break
            end
        end
        plot(1:numel(q),repmat(min(q)+deltaq*a/100,[numel(q),1]),'o')
        areaq=min(q)+deltaq*a/100;
%         figure(5)
%         histogram(q);%,round(sqrt(numel(q))))
        qMeanVec=[qMeanVec;mean(q)];
        qMidVec=[qMidVec;(min(q)+mean(q))/2];
        qBalVec=[qBalVec;areaq];
        figure(2);
        plot(qMeanVec,'r');hold on;
        plot(qMidVec,'g');
        plot(qBalVec,'b');
        disp('Calculating perm')
        K = (Nz/(Nx*Ny))*((meanq*1)/(1-0))/(voxelSize*1)*10^12;
        Kmin=(Nz/(Nx*Ny))*((min(q)*1)/(1-0))/(voxelSize*1)*10^12;
%         Kmidq=(Nz/(Nx*Ny))*(((min(q)+mean(q))/2*1)/(1-0))/(voxelSize*1)*10^12;
        Kbalance=(Nz/(Nx*Ny))*((areaq*1)/(1-0))/(voxelSize*1)*10^12;
        permVec=[permVec;K];
        minPermVec=[minPermVec;Kmin];
%         midqPermVec=[midqPermVec;Kmidq];
        balancePermVec=[balancePermVec;Kbalance];
        disp(['Calculated perm: ',num2str(K)])
        if singlePassFlag
            return 
        end
%         toc
%         tic   
       [dualGridNDFP] = distributeFaceData(numPartitionsX,numPartitionsY,numPartitionsZ,midBCPressures,...
                                         numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,...
                                         dualFaceLocsX,dualFaceLocsY,dualFaceLocsZ);
        %% dual grid routine
        disp('Solving dual fine partitions')
        dualMidBCPressures=cell(numTotDualParts,1);
%         parfor (n=1:numTotDualParts,parForArg)
        for (n=1:numTotDualParts)
            dualTime=tic;
%         for i=1:numDualPartitionsX
%             for j=1:numDualPartitionsY
%                 for k=1:numDualPartitionsZ
            [i,j,k]=ind2sub([numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ],n);
                    if rmDualFlags(i,j,k)==0
                        solver='AGMG';
                    else
                        solver='MLDIVIDE';
                    end
                    localOrigin=[dualFaceLocsX(i,1),dualFaceLocsY(j,1),dualFaceLocsZ(k,1)];
                    if ~diskFlag
                        S=transDualParts{i,j,k};
                        if ~microFlag
                            G3=Gdualparts{i,j,k};
                        else
                            dimVec=[i==1||i==numDualPartitionsX,...
                                    j==1||j==numDualPartitionsY,...
                                    k==1||k==numDualPartitionsZ];
                            dimMask=[numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ];
                            dimMask=dimMask>1;
                            dimVec=dimVec.*dimMask;
                            if sum(dimVec)==3
                                G3=G8;
                            elseif sum(dimVec)==0
                                G3=G1; 
                            elseif sum(dimVec)==1
                                [~,gridInd]=max(dimVec);
                                G3=GH{gridInd};
                            elseif sum(dimVec)==2
                                [~,gridInd]=min(dimVec);
                                G3=G4{gridInd}; 
                            end
                        end
                    else
                        disp(['Loading dual grid ',num2str(i),num2str(j),num2str(k),' to RAM'])
                        S=serialisedloadsave([diskDir,'dualTrans',num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        G3=serialisedloadsave([diskDir,'dualGrid',num2str(i),num2str(j),num2str(k),'.mat'],[],0,diskDumpSerialiserFlag);
                        if ~diskDumpSerialiserFlag
                            S=S.data;
                            G3=G3.data;
                        end
                    end
                    % determine active faces
                    activeFaces=[i==1,i==numDualPartitionsX;
                                 j==1,j==numDualPartitionsY;
                                 k==1,k==numDualPartitionsZ];
                    if sum(activeFaces(:,1))==3 ||sum(activeFaces(:,2))==3
                        dualFlag=true; % corners do not contribute
                    else
                        dualFlag=false;
                    end
                    % simply the centre of the current grid
                    localDual=[mean(dualFaceLocsX(i,:))-localOrigin(1) mean(dualFaceLocsX(i,:))-localOrigin(1);
                               mean(dualFaceLocsY(j,:))-localOrigin(2) mean(dualFaceLocsY(j,:))-localOrigin(2);
                               mean(dualFaceLocsZ(k,:))-localOrigin(3) mean(dualFaceLocsZ(k,:))-localOrigin(3)];
                    % NDFP
                    NDFP=dualGridNDFP{i,j,k};
                    dualBc=[dualFaceLocsX(i,:);
                            dualFaceLocsY(j,:);
                            dualFaceLocsZ(k,:)]-localOrigin';
                    [~,Pressure] =nonDirichletCorrectionSubroutine3D...
                        (G3,S,NDFP,dualBc,localOrigin,localDual,plotFlag,iterNum,...
                        dualFlag,rmDualFlags(i,j,k),~activeFaces,[],[]);
                     dualMidBCPressures{n}=Pressure;
                      disp(['Solved dual fine partition ',num2str(i),num2str(j),num2str(k),...
                          ' with ', solver,' ',num2str(toc(dualTime)),': seconds elapsed'])
%                 end
%             end
%         end
        end
        clear G3
%         toc
%         tic
       [gridNDFP] = distributeFaceData(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,dualMidBCPressures,...
                                         numPartitionsX,numPartitionsY,numPartitionsZ,...
                                         faceLocsX,faceLocsY,faceLocsZ);
        disp('Calculating residual')
%         residual=max(abs(diff(q)))/mean(q);
        residual=(abs(max(q))-abs(min(q)))/mean(q);
        %%%%%%%%%%%%%%
        resVec=[resVec;residual];
        timeVec=[timeVec;toc(simTime)];
        figure(3)
        semilogy(timeVec,resVec);
        drawnow
        figure(4)
        plot(timeVec,permVec);hold on;
        plot(timeVec,minPermVec);
        midPermVec=(minPermVec+permVec)/2;
        plot(timeVec,midPermVec);
        figure(5)
        plot(timeVec,permVec);hold on;
        plot(timeVec,minPermVec);
        midPermVec=(minPermVec+permVec)/2;
        plot(timeVec,midPermVec);
        plot(timeVec,balancePermVec,'o')
%         plot(timeVec,midqPermVec,'.')
        yi=nan;
        yj=nan;
        yk=nan;
        if iterNum>1
            [x1,y1] = plotLine(timeVec(end-1:end),permVec(end-1:end),[timeVec(end-1) timeVec(end-1)*4]);
            [x2,y2] = plotLine(timeVec(end-1:end),minPermVec(end-1:end),[timeVec(end-1) timeVec(end-1)*4]);
            [x3,y3] = plotLine(timeVec(end-1:end),midPermVec(end-1:end),[timeVec(end-1) timeVec(end-1)*4]);
            [xi,yi] = polyxpoly(x1,y1,x2,y2);
            [xj,yj] = polyxpoly(x1,y1,x3,y3);
            [xk,yk] = polyxpoly(x3,y3,x2,y2);
            hold on
            plot([xi,xj,xk],[yi,yj,yk],'*');
            hold off
        end
        
        title({['Midpoint K: ',num2str((minPermVec(end)+permVec(end))/2)],...
               ['Extrap K: ',num2str((yi+yk+yk)/3)],...
               ['Balance K: ',num2str(balancePermVec(end))]})
%         [asympK] = extrapolateAsymptoteYDW(permVec);
        %         figure(5)
%         plot(1:numel(timeVec),timeVec)
        drawnow
        disp(['Calculated perm: ',num2str(K)])
%         disp(['Extrapolated perm: ',num2str(asympK)])

        disp(['Calculated residual: ', num2str(residual)])
        disp(['Dual cycle iteration number: ', num2str(iterNum)])
        disp(['Time elapsed: ', num2str(toc)])
        if plotFlag
%             clf(1)   
        end
        fprintf(fileID,'%d %d %d %d\r\n',[iterNum toc K residual]);
        fclose(fileID);
    end
%     save(['vels',num2str(iterNum),'.mat'],'velDat','-v7.3')
%     profile report
%     profile off
data.permVec=permVec;
data.timeVec=timeVec;
data.resVec=resVec;
data.minPermVec=minPermVec;
end

function [dualGridNDFP] = distributeFaceData(numPartitionsX,numPartitionsY,numPartitionsZ,midBCPressures,...
                                         numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,...
                                         dualFaceLocsX,dualFaceLocsY,dualFaceLocsZ)
       disp('Aligning pressures') % embed this algorithm into the loop later if you want?
        leftFaceData=[];%zeros(Ny*Nz*numPartitionsX+Nx*Nz*numPartitionsY+Nx*Ny*numPartitionsZ,1);
        rightFaceData=[];%zeros(Ny*Nz*numPartitionsX+Nx*Nz*numPartitionsY+Nx*Ny*numPartitionsZ,1);
%         centralFaceData=[];
        for i=1:numPartitionsX
            for j=1:numPartitionsY
                for k=1:numPartitionsZ
                    n=sub2ind([numPartitionsX,numPartitionsY,numPartitionsZ],i,j,k);
                    if ~isempty(midBCPressures{n})
                        leftFaceData=[leftFaceData;
                                       midBCPressures{n}.left.centroids midBCPressures{n}.left.pressure];
                        rightFaceData=[rightFaceData;
                                       midBCPressures{n}.right.centroids midBCPressures{n}.right.pressure];     
    %                     centralFaceData=[centralFaceData;
    %                                    midBCPressures{i,j,k}.centroids midBCPressures{i,j,k}.pressure];
                    end
                end
            end
        end
        dualGridNDFP=cell(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ);
        for i=1:numDualPartitionsX
            for j=1:numDualPartitionsY
                for k=1:numDualPartitionsZ
                    boundingBox=[dualFaceLocsX(i,:);
                                 dualFaceLocsY(j,:);
                                 dualFaceLocsZ(k,:)];
                    %minimums need RHS, maximums need LHS, active faces
                    %take dummy data
                    inactiveFaces=[i==1,i==numDualPartitionsX;
                                   j==1,j==numDualPartitionsY;
                                   k==1,k==numDualPartitionsZ];
                    NDFP=cell(3,2);
                    for m=1:3
                        for n=1:2
                            otherDims=setdiff([1,2,3],m);
%                             otherBound=setdiff([1,2],n);
                            if inactiveFaces(m,n)==1
                                NDFP{m,n}=[];
                                continue
                            end
                            if n==1
                                tempInds=rightFaceData(:,m)==boundingBox(m,n) &...% squeeze along m dim
                                         rightFaceData(:,otherDims(1))>boundingBox(otherDims(1),1) &...; % bound along other dims
                                         rightFaceData(:,otherDims(1))<boundingBox(otherDims(1),2) &...
                                         rightFaceData(:,otherDims(2))>boundingBox(otherDims(2),1) &...
                                         rightFaceData(:,otherDims(2))<boundingBox(otherDims(2),2);
                                NDFP{m,n}=rightFaceData(tempInds,:);
                            elseif n==2
                                tempInds=leftFaceData(:,m)==boundingBox(m,n) &...% squeeze along m dim
                                         leftFaceData(:,otherDims(1))>boundingBox(otherDims(1),1) &...; % bound along other dims
                                         leftFaceData(:,otherDims(1))<boundingBox(otherDims(1),2) &...
                                         leftFaceData(:,otherDims(2))>boundingBox(otherDims(2),1) &...
                                         leftFaceData(:,otherDims(2))<boundingBox(otherDims(2),2);
                                NDFP{m,n}=leftFaceData(tempInds,:);
                            end
                        end
                    end
                    dualGridNDFP{i,j,k}=NDFP;
                end
            end
        end     
end

function [x,y] = plotLine(A,B,xlim)
    m = (B(2)-B(1))/(A(2)-A(1));
    n = B(2) - A(2)*m;
    y1 = m*xlim(1) + n;
    y2 = m*xlim(2) + n;
    hold on
    line([xlim(1) xlim(2)],[y1 y2])
    plot(A,B,'*')
    x=[xlim(1) xlim(2)];
    y=[y1 y2];
    hold off
end
%%
