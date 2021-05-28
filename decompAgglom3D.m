%% 3d version of domain decomp - unvectorised parameter extension - no onverlapping dual grid
%% things to improve:
% velocity fields should be calculated in the loop for partition flow rate,
% not recovered outside the loop.
% 

function varargout = decompAgglom3D(domainFile,alpha,voxelSize,numPartitions,agglomFlag,preProcDiskDumpFlag,simulationDiskDumpFlag,runPartsInSerial)
    %% Agglomeration Module
    %% extra data for 3D
    data=struct();
    resTol=1e-2;
    numPartitionsX=numPartitions(1);
    numPartitionsY=numPartitions(2);
    numPartitionsZ=numPartitions(3);
    incX=double(numPartitionsX>1);
    incY=double(numPartitionsY>1);
    incZ=double(numPartitionsZ>1);
    singlePassFlag=numPartitionsX+numPartitionsY+numPartitionsZ==3;    
    %% load binary
    %%
    disp('Loading Data')
    domain=load([pwd,'/',domainFile,'.mat'],'-mat');
    domain=domain.(domainFile);
    domain=domain(1:128,1:128,1:128);
    [Nx, Ny, Nz] = size(domain);
    dualLengthZ=ceil(Nz/numPartitionsZ/2);
    dualLengthY=ceil(Ny/numPartitionsY/2); % allows for odd numbers. facility to safely overlap is not implemented
    dualLengthX=ceil(Nx/numPartitionsX/2);
        
    %% prepare interim data file
    fileID = fopen('OutputLog.txt','a');
    fprintf(fileID,['\r\n\r\n Decompositional Algorithm 3D running for new test on ',domainFile,...
                    ' sized ',num2str(Nx),' ',num2str(Ny),' ',num2str(Nz) ,...
                    ' split into ',num2str(numPartitions),' parts with dual visibility of ',num2str(dualLengthZ), ' ', num2str(dualLengthY), ' ', num2str(dualLengthZ),' on ',datestr(datetime),'\r\n\r\n']);
    %% set secondary flags
    if preProcDiskDumpFlag==true
        disp('Preprocessing Disk Dumping is Active')
    end
    rootDir=pwd;
    diskDir=[rootDir,'/transDump/'];
    delta=0;
    if runPartsInSerial
      parForArg = 0;
    else
      parForArg = Inf;
    end
    numDualPartitionsZ=numPartitionsZ+incZ;
    numDualPartitionsY=numPartitionsY+incY;
    numDualPartitionsX=numPartitionsX+incX;
    %% remove cells
    [domain,~] = removeDisconnections(domain);
    %% calculate weights
    disp('Calculating fdgpa conductivity')
    D=bwdist(domain);
%     alpha=1;
    maxdist=max(D(:));
    if agglomFlag
        bwD=D;
    end
    [Dmax] = FindWW2(D, maxdist,Nx, Ny, Nz, alpha,voxelSize);
    D=alpha.*voxelSize.^2./8.*(2.*Dmax.*D-D.^2);
    if agglomFlag
        D=D.*voxelSize;
    else
        clear Dmax
    end
    %% anchor image with eps pseudo connectivity
    anchors=false(size(D));
    domain=D==0;
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
    if ~singlePassFlag
        disp('Seeding domain with anchors')       
        for i=1:numPartitionsX
            for j=1:numPartitionsY
                for k=1:numPartitionsZ
                    temp=domain(intervalX{i},intervalY{j},intervalZ{k});
                    localOrigin=[min(intervalX{i})-1,min(intervalY{j})-1,min(intervalZ{k})-1];
                    %disp(['Seeding Tendrils for grid block',num2str(i),num2str(j),num2str(k)])
                    [anchorPoints] = seedTendrils(temp);
                    if ~isempty(anchorPoints)
                        for m=1:size(anchorPoints,1)
                            x=anchorPoints(m,1)+localOrigin(1);
                            y=anchorPoints(m,2)+localOrigin(2);
                            anchors(x,y,:)=true;
                        end
                    end
                end
            end
        end
        for i=1:numDualPartitionsX
            for j=1:numDualPartitionsY
                for k=1:numDualPartitionsZ
                    temp=domain(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
                    localOrigin=[min(dualIntervalX{i})-1,min(dualIntervalY{j})-1,min(dualIntervalZ{k})-1];
                    %disp(['Seeding Tendrils for dual grid block',num2str(i),num2str(j),num2str(k)])
                    [anchorPoints] = seedTendrils(temp);
                    if ~isempty(anchorPoints)
                        for m=1:size(anchorPoints,1)
                            x=anchorPoints(m,1)+localOrigin(1);
                            y=anchorPoints(m,2)+localOrigin(2);
                            anchors(x,y,:)=true;
                        end
                    end
                end
            end
        end       
    %     disp('Superimposing anchors with weights')
        D(anchors)=D(anchors)+eps;
    end
    clear temp anchors domain
    %% split weights and binary and NDCF arrays
%     disp('Splitting coarse partitions')
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
    %%
    gridNDFP=cell(numPartitionsX,numPartitionsY,numPartitionsZ);
    rmFlags=zeros(numPartitionsX,numPartitionsY,numPartitionsZ);
    for i=1:numPartitionsX
        for j=1:numPartitionsY
            for k=1:numPartitionsZ
                rock=struct();
                disp(['Splitting primal partition ', num2str(i),num2str(j),num2str(k)])
                %disp('Generating representative coarse grid')
                temp=D(intervalX{i},intervalY{j},intervalZ{k});
                    G=[];
                    %disp('Microporosity inactive, generating reduced cartesian grid set for each partition')
                    [~,rmFlag] = removeDisconnections(~temp);
                    if rmFlag
                        error('Disconnected subdomain detected.')
                    end
                    rmFlags(i,j,k)=rmFlag;
                    if agglomFlag
%                         tempDm=DDown2(intervalX{i}(2:2:end)./2,intervalY{j}(2:2:end)./2,intervalZ{k}(2:2:end)./2);
                        bwTemp=bwD(intervalX{i},intervalY{j},intervalZ{k});
                        tempDm=Dmax(intervalX{i},intervalY{j},intervalZ{k});
%                         locUpMap=UpMap(intervalX{i}(2:2:end)./2,intervalY{j}(2:2:end)./2,intervalZ{k}(2:2:end)./2);
                        localOrigin=[faceLocsX(i,1),faceLocsY(j,1),faceLocsZ(k,1)];
                        [Gparts{i,j,k}, rock] = generateAgglomeratedGridYDW(bwTemp,temp,tempDm,localOrigin, voxelSize);
                    else
                        rock.perm=double(temp(:))+delta;
                        indexMap=(1:numel(rock.perm))';
                        Gparts{i,j,k}=generateReducedGridYDW(size(temp),indexMap(rock.perm==0|isnan(rock.perm)));
                        rock.perm=rock.perm(rock.perm>0);
                        Gparts{i,j,k}=mcomputeGeometry(Gparts{i,j,k});
                    end
                    transParts{i,j,k}=computeTrans(Gparts{i,j,k},rock);
                    if preProcDiskDumpFlag
                        disp('Disk dumping active, dumping trans and grid data to disk')
                        serialisedloadsave([diskDir,'trans', num2str(i),num2str(j),num2str(k),'.mat'],transParts{i,j,k},1,1);
                        transParts{i,j,k}=[];
                        serialisedloadsave([diskDir,'grid', num2str(i),num2str(j),num2str(k),'.mat'],Gparts{i,j,k},1,1);
                        Gparts{i,j,k}=[];
                    end
            end
        end
    end
    %% split into dual arrays
%     disp('Splitting dual coarse partitions')
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
                    disp(['Splitting dual partition ', num2str(i),num2str(j),num2str(k)])
                    %disp('Generating representative dual coarse grid')
                    temp=D(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
%                         disp('Microporosity inactive, generating reduced cartesian grid set for each dual partition')
                        [~,rmDualFlag] = removeDisconnections(~temp);
                        if rmDualFlag
                            error('Disconnected subdomain detected.')
                        end
                        rmDualFlags(i,j,k)=rmDualFlag;
                        if agglomFlag
%                             tempDm=DDown2(dualIntervalX{i}(2:2:end)./2,dualIntervalY{j}(2:2:end)./2,dualIntervalZ{k}(2:2:end)./2);
                            bwTemp=bwD(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
                            tempDm=Dmax(dualIntervalX{i},dualIntervalY{j},dualIntervalZ{k});
%                             locUpMap=UpMap(dualIntervalX{i}(2:2:end)./2,dualIntervalY{j}(2:2:end)./2,dualIntervalZ{k}(2:2:end)./2);
                            localOrigin=[dualFaceLocsX(i,1),dualFaceLocsY(j,1),dualFaceLocsZ(k,1)];
                            [Gdualparts{i,j,k}, rock] = generateAgglomeratedGridYDW(bwTemp,temp,tempDm,localOrigin,voxelSize);
                        else
                            rock.perm=double(temp(:))+delta;
                            indexMap=(1:numel(rock.perm))';
                            Gdualparts{i,j,k}=generateReducedGridYDW(size(temp),indexMap(rock.perm==0|isnan(rock.perm)));
                            rock.perm=rock.perm(rock.perm>0);
                            Gdualparts{i,j,k}=mcomputeGeometry(Gdualparts{i,j,k});
                        end
                        transDualParts{i,j,k}=computeTrans(Gdualparts{i,j,k},rock); 
                        if preProcDiskDumpFlag
                            disp('Disk dumping active, dumping trans and grid data to disk')
                            serialisedloadsave([diskDir,'dualTrans', num2str(i),num2str(j),num2str(k),'.mat'],transDualParts{i,j,k},1,1);
                            transDualParts{i,j,k}=[];
                            serialisedloadsave([diskDir,'dualGrid', num2str(i),num2str(j),num2str(k),'.mat'],Gdualparts{i,j,k},1,1);
                            Gdualparts{i,j,k}=[];
                        end
                end
            end
        end
    end
    %% face locations
%     disp('Setting face locations') % regroup boundary face locations to suit grids
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
    %% reload grid and trans data if flags
    if preProcDiskDumpFlag==true
        if simulationDiskDumpFlag==false
            disp('Loading all partition data back to RAM')
            for i=1:numPartitionsX
                for j=1:numPartitionsY
                    for k=1:numPartitionsZ
                        Gparts{i,j,k}=serialisedloadsave([diskDir,'grid', num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        transParts{i,j,k}=serialisedloadsave([diskDir,'trans', num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        if ~1
                            Gparts{i,j,k}=Gparts{i,j,k}.data;
                            transParts{i,j,k}=transParts{i,j,k}.data;
                        end
                        disp(['Primal partition ', num2str(i),num2str(j),num2str(k),' loaded back into RAM'])
                    end
                end
            end
            for i=1:numDualPartitionsX
                for j=1:numDualPartitionsY
                    for k=1:numDualPartitionsZ
                        Gdualparts{i,j,k}=serialisedloadsave([diskDir,'dualGrid',num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        transDualParts{i,j,k}=serialisedloadsave([diskDir,'dualTrans',num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        if ~1
                            Gdualparts{i,j,k}=Gdualparts{i,j,k}.data;
                            transDualParts{i,j,k}=transDualParts{i,j,k}.data;
                        end
                        disp(['Dual partition ',num2str(i),num2str(j),num2str(k),' loaded back into RAM'])
                    end
                end
            end
        end
    end
    %%  solver routine
    % iterate until matching!
    % generate coarse grids
    disp('Solving...')
    residual=100;
    iterNum=0;
    resVec=[];
    permVec=[];
    minPermVec=[];
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
    numIterations = 50;
    stateMat=cell(numPartitionsX,numPartitionsY,numPartitionsZ);

    while iterNum<=numIterations && residual>resTol
        fileID = fopen('OutputLog.txt','a');
        midBCPressures=cell(numTotParts,1);
        subVel=cell(numTotParts,1);
        subVelZ=cell(numTotParts,1);
        subVelX=cell(numTotParts,1);
        subVelY=cell(numTotParts,1);
        sumVelZ=cell(numTotParts,1);
        %solve flow with coarse grids
        iterNum=iterNum+1;
        disp('Solving primal partitions')
%         toc
%         tic
        localDual=[cellDualFaceLocsX(1,:);
                   cellDualFaceLocsY(1,:);
                   cellDualFaceLocsZ(1,:)];
        clear totVel

        %% grid routine
%         parfor (n=1:numTotParts,parForArg)
        for (n=1:numTotParts)
            partTime=tic;
            [i,j,k]=ind2sub([numPartitionsX,numPartitionsY,numPartitionsZ],n);
                    if rmFlags(i,j,k)==0
                        solver='AGMG';
                    else
                        solver='MLDIVIDE';
                    end
                    localOrigin=[faceLocsX(i,1),faceLocsY(j,1),faceLocsZ(k,1)];
                    dualFlag=false;
                    if ~simulationDiskDumpFlag
                        S=transParts{i,j,k};
                            G2=Gparts{i,j,k};
                    else
                        disp(['Loading primal grid ',num2str(i),num2str(j),num2str(k),' to RAM'])
                        S=serialisedloadsave([diskDir,'trans',num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        G2=serialisedloadsave([diskDir,'grid',num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        if ~1
                            G2=G2.data;
                            S=S.data;
                        end
                    end
                    if iterNum==1 % populate NDFP with initial guess data
                            NDFP{1,1}=dirPressVecX(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{1,2}=dirPressVecX(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{2,1}=dirPressVecY(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{2,2}=dirPressVecY(round(localOrigin(3)+(localDual(3,1)+localDual(3,2))/2));
                            NDFP{3,1}=dirPressVecZ(k);
                            NDFP{3,2}=dirPressVecZ(k+1); 
                        init=[];
                    else
                        NDFP=gridNDFP{i,j,k};
                        init=[];
                    end
                    % determine active faces
                    activeFaces=[i==1,i==numPartitionsX;
                                 j==1,j==numPartitionsY;
                                 k==1,k==numPartitionsZ];
                    [state,dualPressure] =nonDirichletCorrectionSubroutine3D...
                        (G2,S,NDFP,bc,localOrigin,localDual,1,iterNum,...
                        dualFlag,rmFlags(i,j,k),~activeFaces,init,[],agglomFlag);
                    midBCPressures{n}=dualPressure;
                    disp(['Solved primal partition ',num2str(i),num2str(j),num2str(k),...
                          ' with ', solver,' ',num2str(toc(partTime)),': seconds elapsed'])
        %             disp('Extracting Flux')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ADD option to collect velP fields%%%%
                    if nargout==2 && ~singlePassFlag
                        disp('Velocity and Pressure Fields Requested: Collating.')
                        stateMat{i,j,k}=state;
                    end
                    tempVel=faceFlux2cellVelocity(G2,state.flux);
                    tempVel=tempVel(:,3);
                    if agglomFlag % getting the agglomerated q cross section
                        velData = [G2.cells.centroids, tempVel];
%                         [zCoords, indsa, indsb] = unique(velData(:,3));
%                         zVels=accumarray(indsb, velData(indsb));
                        rc = unique(velData(:,3));
                        zCoords = zeros(size(rc,1),2);
                        for x = 1:size(rc,1)
                            xxx = rc(x);
                            temp = velData(:,3) == xxx;
                            temp2 =  sum(velData(temp,4));
                            zCoords(x,:)= [round(xxx,1), temp2 ];
                        end
                        t = 0;
                        for m= 1:2:max(zCoords(:,1))+0.5
                            t = t + 1;
                            q1 = round(zCoords(:,1),1) == m;
                            q11 = (zCoords(q1,2));
                            q2 = round(zCoords(:,1),1) == m-0.5;
                            q12 = zCoords(q2,2);
                            q3 = round(zCoords(:,1),1) == m+0.5;
                            q13 = zCoords(q3,2);
                            if  isempty(q11)
                                q11 = 0;
                            end
                            if  isempty(q12)
                                q12 = 0;
                            end
                            if  isempty(q13)
                                q13 = 0;
                            end
                            q14(t) = q11+q12+q13;
                        end
                        %clear q1 q2 q3 q11 q12 q13 face1 faceIndexList
                        q14 = q14(:);
                        q14=permute(q14,[3,2,1]);
                        sumVelZ{n}=q14;
                    else
                        cartVelZ=zeros(G2.cartDims);
                        cartVelZ=cartVelZ(:);
                        cartVelZ(G2.cells.indexMap)=tempVel;
                        tempVel=reshape(cartVelZ,G2.cartDims);
                        tempuz=zeros(1,1,G2.cartDims(3));
                        for (p = 1:G2.cartDims(3))
                            tempuz(1,1,p) = sum(sum(tempVel(:,:,p)));
                        end
                        sumVelZ{n}=tempuz;
                    end
        end
        %%
        clear G2
        disp('Reconstructing Z Fluxes')
        sumVelZ=reshape(sumVelZ,[numPartitionsX,numPartitionsY,numPartitionsZ]);
        sumVelZ=cell2mat(sumVelZ);
        tempuz = zeros(size(sumVelZ,3),1);
        for (k = 1:size(sumVelZ,3))
            tempuz(k) = sum(sum(sumVelZ(:,:,k)));
        end
        if agglomFlag
            q = tempuz(:)./2;
        else
            q = tempuz(:).*voxelSize;
        end
        q=q';
        meanq = mean(q);
        figure(1);clf(1);
        title('Flow Z')
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
            varargout{1}=balancePermVec(end);
            if nargout == 2
                % collate into cart grids!
                varargout{2}=state;
            end
            return 
        end
%         toc
%         tic   
        if ~agglomFlag
            [dualGridNDFP] = distributeFaceData(numPartitionsX,numPartitionsY,numPartitionsZ,midBCPressures,...
                                         numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,...
                                         dualFaceLocsX,dualFaceLocsY,dualFaceLocsZ);
        else
            [dualGridNDFP] = distributeAgglomeratedFaceData(numPartitionsX,numPartitionsY,numPartitionsZ,midBCPressures,...
                                         numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,...
                                         dualFaceLocsX,dualFaceLocsY,dualFaceLocsZ);
        end
        %% dual grid routine
        disp('Solving dual partitions')
        dualMidBCPressures=cell(numTotDualParts,1);
        parfor (n=1:numTotDualParts,parForArg)
%         for (n=1:numTotDualParts)
            dualTime=tic;
            [i,j,k]=ind2sub([numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ],n);
                    if rmDualFlags(i,j,k)==0
                        solver='AGMG';
                    else
                        solver='MLDIVIDE';
                    end
                    localOrigin=[dualFaceLocsX(i,1),dualFaceLocsY(j,1),dualFaceLocsZ(k,1)];
                    if ~simulationDiskDumpFlag
                        S=transDualParts{i,j,k};
                            G3=Gdualparts{i,j,k};
                    else
                        disp(['Loading dual grid ',num2str(i),num2str(j),num2str(k),' to RAM'])
                        S=serialisedloadsave([diskDir,'dualTrans',num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
                        G3=serialisedloadsave([diskDir,'dualGrid',num2str(i),num2str(j),num2str(k),'.mat'],[],0,1);
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
                        (G3,S,NDFP,dualBc,localOrigin,localDual,1,iterNum,...
                        dualFlag,rmDualFlags(i,j,k),~activeFaces,[],[],agglomFlag);
                     dualMidBCPressures{n}=Pressure;
                      disp(['Solved dual partition ',num2str(i),num2str(j),num2str(k),...
                          ' with ', solver,' ',num2str(toc(dualTime)),': seconds elapsed'])
%                 end
%             end
%         end
        end
        %%
        clear G3
%         toc
%         tic
        if ~agglomFlag
            [gridNDFP] = distributeFaceData(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,dualMidBCPressures,...
                                         numPartitionsX,numPartitionsY,numPartitionsZ,...
                                         faceLocsX,faceLocsY,faceLocsZ);
        else
            [gridNDFP] = distributeAgglomeratedFaceData(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,dualMidBCPressures,...
                                         numPartitionsX,numPartitionsY,numPartitionsZ,...
                                         faceLocsX,faceLocsY,faceLocsZ);
        end
%         disp('Calculating residual')
%         residual=max(abs(diff(q)))/mean(q);
        residual=(abs(max(q))-abs(min(q)))/mean(q);
        %%%%%%%%%%%%%%
        resVec=[resVec;residual];
        timeVec=[timeVec;toc(simTime)];
        figure(2)
        title('Perm vs Time')
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
        disp(['Time elapsed: ', num2str(timeVec(end))])

        fprintf(fileID,'%d %d %d %d\r\n',[iterNum timeVec(end) K residual]);
        fclose(fileID);
    end

if nargout==2
    fullState=struct();
    disp('Velocity and Pressure Fields Requested: Final Domain Stitching.')
    for i=1:numPartitionsX
        for j=1:numPartitionsY
            for k=numPartitionsZ
                temp=stateMat{i,j,k};
                % need G as well
            end
        end
    end
end
varargout{1}=balancePermVec(end);
varargout{2}=fullState;
end

function [dualGridNDFP] = distributeAgglomeratedFaceData(numPartitionsX,numPartitionsY,numPartitionsZ,midBCPressures,...
                                         numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,...
                                         dualFaceLocsX,dualFaceLocsY,dualFaceLocsZ)
%     disp('Aligning pressures') % embed this algorithm into the loop later if you want?
    midBCPressures=cell2mat(midBCPressures); %cell2mat is cancer, avoid it
    midBCPressures=unique(midBCPressures,'rows');
    %for each dual grid, find the boundaries
    dualGridNDFP=cell(numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ);
    for i=1:numDualPartitionsX
        for j=1:numDualPartitionsY
            for k=1:numDualPartitionsZ
                primalFaceCentroids=midBCPressures(:,[1:3]); 
                primalFacePressures=midBCPressures(:,4); 
                % for each of the 6 faces, find the upscaled face
                % values in the order [xmin xmax; ymin ymax; zmin zmax]
                boundingBox=[dualFaceLocsX(i,:);
                             dualFaceLocsY(j,:);
                             dualFaceLocsZ(k,:)];

                inactiveFaces=[i==1,i==numDualPartitionsX;
                               j==1,j==numDualPartitionsY;
                               k==1,k==numDualPartitionsZ];
                % extract the subset from the upscaleFaceCentroids
                locFaceCentroids = primalFaceCentroids(primalFaceCentroids(:,1)>=boundingBox(1,1) & ...
                                             primalFaceCentroids(:,1)<=boundingBox(1,2) & ...
                                             primalFaceCentroids(:,2)>=boundingBox(2,1) & ...
                                             primalFaceCentroids(:,2)<=boundingBox(2,2) & ...
                                             primalFaceCentroids(:,3)>=boundingBox(3,1) & ...
                                             primalFaceCentroids(:,3)<=boundingBox(3,2),:);
                locFacePressures = primalFacePressures(primalFaceCentroids(:,1)>=boundingBox(1,1) & ...
                                             primalFaceCentroids(:,1)<=boundingBox(1,2) & ...
                                             primalFaceCentroids(:,2)>=boundingBox(2,1) & ...
                                             primalFaceCentroids(:,2)<=boundingBox(2,2) & ...
                                             primalFaceCentroids(:,3)>=boundingBox(3,1) & ...
                                             primalFaceCentroids(:,3)<=boundingBox(3,2),:);

                NDFP=cell(3,2);
                for m=1:3 %find primal faces, then map each face to the upscaled value [no left and right nonsense]
                    for n=1:2
                        otherDims=setdiff([1,2,3],m);
%                             otherBound=setdiff([1,2],n);
                        if inactiveFaces(m,n)==1
                            NDFP{m,n}=[];
                            continue
                        end
                        %find the faces on the boundary
                        locFaceInds=locFaceCentroids(:,m)==boundingBox(m,n); %retains ordering built in G
                        tempFaceCentroids=locFaceCentroids(locFaceInds,:); 
                        tempFacePressures=locFacePressures(locFaceInds);
                        NDFP{m,n}=[tempFaceCentroids, tempFacePressures];
                        NDFP{m,n}=unique(NDFP{m,n},'stable', 'rows');
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
                dualGridNDFP{i,j,k}=NDFP;
            end
        end
    end                        
end

function [dualGridNDFP] = distributeFaceData(numPartitionsX,numPartitionsY,numPartitionsZ,midBCPressures,...
                                         numDualPartitionsX,numDualPartitionsY,numDualPartitionsZ,...
                                         dualFaceLocsX,dualFaceLocsY,dualFaceLocsZ)
%     disp('Aligning pressures') % embed this algorithm into the loop later if you want?
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
