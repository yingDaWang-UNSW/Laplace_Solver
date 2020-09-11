
% %% startup mrst2016bydw
% addpath(genpath(pwd))
% startup
% gravity reset on
% mrstModule add libgeometry incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
% mrstModule add coarsegrid
% mrstModule add agmg
% addpath(genpath(pwd))
% % cd ..
%% generate datalist and solutionlist
% dataList={'geopack';'bentheimer';'Berea';'carbbig2';'geoclear';'geomin';'ROIsandporesser'};
dataList={'yulai100', 'yulai200', 'yulai300', 'yulai400', 'yulai600', 'yulai1200'};

solnList=cell(numel(dataList),1);
%% generate input arguments
% voxelSizeList=[10;5.3;5.3;3.8;16.5;16.5].*1e-6;
% alphaList=[1;1;1;2;2;1];
voxelSizeList=[45.6, 22.8, 15.2, 11.4, 7.6, 3.8].*1e-6;
alphaList=[1, 1, 1, 1, 1, 1];
%trunc: how much of the input image you would like to use in the 3 basis directions
%voxelSize: the true size of each voxel in m
%alpha: the fracture or pore shape factor
%YDFlag: whether to use the YD or TC Dmax algorithm
%numPartitions: the number of subdomains [must be a factor of truncz]
%dualLength: the size of the dual grid [must be between 1 and truncz/numpartitions-1] [non-overlap is at truncz/numpartitions/2]
%relaxFlag: a flag that determines if the outer iteration applies a relaxation factor to the sequential input [overrelaxation can improve convergence, but there exists an upper stability limit]
%precondFlag: a flag that determines if the inner iteration uses the previous results of the outer iteration as the initial estimate.
%pressureBasisFlag: a flag that determines if the first iteration is performed with a boundary approximation obtained from an upscaled solution to the overall domain.
%microFlag: a flag that determines if the dmoain is solved in its entirety,
%or if the solid elements are completely impermeable. 
%agmgFlag: a flag that determines if the agmg solver is used to solve the
%system. If this is false, then no anchoring will occur and the system will
%solve using a mixture of MLDIVIDE and AGMG depending on natural anchoring
%unixflag: self explanatory
%diskFlag: whether the algorithm will dump and retrieve grid and perm data from disk to save RAM
%runInSerial: whether the solver will solve partitions in parallel or series
%plotFlag: will plot the solution of flux and pressure in real time within
%the inner loop. Best visual results if agmg is deactivated. Only works if
%serialised
%% run loop note: in this version, the input arguments are varied algorithmically in the loop itself. you may want to pre define them in arrays.
numTruncVec=[100, 100, 100;
             200, 200, 200;
             300, 300, 300;
             400, 400, 400;
             600, 600, 600;
             1200, 1200, 1200;];
            
numPartitionsVec=[1, 1, 1;
                  1, 1, 1;
                  1, 1, 1;
                  1, 1, 1;
                  1, 1, 1;
                  4, 4, 4]; % the domain must be divisible by double the partition number, or 4x if agglomeration is active.

for i=5:numel(solnList)
%     truncx=1:50;%size(image,1);
%     truncy=1:50;%size(image,2);
%     %%%%%NOTE ensure truncz length is divisible by numpartitions*2
%     truncz=1:50;%size(image,3);
 
%     for j=1:size(numPartitionsVec,1)
    j=i;
        truncx=1:numTruncVec(i,1);%1210;%size(image,1);%numTruncVec(1);%
        truncy=1:numTruncVec(i,2);%1210;%size(image,2);
        truncz=1:numTruncVec(i,3);%1320;%size(image,3);
        dCalcFlag=true;
        dFile=[dataList{i},'3d'];
        dSaveFlag=false;
        numPartitions=numPartitionsVec(j,:);%20;
        dualLengthCentre = ceil(max(truncz(end)-truncz(1)+1)/numPartitions(3)/2);
        dualLengthOffset=0;%dualLength-1;%number can be between -(dualLength-2) to dualLength-1 
        dualLengthOffsetVec=0;%-(dualLengthCentre-2):5:dualLengthCentre-1;
        for k=1:numel(dualLengthOffsetVec)
            dualLength=dualLengthCentre+dualLengthOffsetVec(k);
            numIterations=50;
            relaxFlag=false;
            relaxIndex=3;
            precondFlag=false;
            pressureBasisFlag=false;
            upscalingFactor=2;
            agglomFlag=false;
            microFlag=false;
            agmgFlag=true; %most issues with agmg have been solved. only deactivate this if plot flag is on and you need cross section face data
            unixFlag=true;         
            condDumpSerialiserFlag=false;
            diskDumpSerialiserFlag=false;
            
            preProcDiskDumpFlag=false;
            simulationDiskDumpFlag=false;
            
            YDFlag=false;
            runYDInSerial=true;
            runPartsInSerial=true;

    %         delete(gcp('nocreate'))
    %         parpool('local',16)
            plotFlag=false;
        %     clf(1);clf(2);clf(3);clf(4);
%             figure(3)
%             hold on
%             figure(4)
%             hold on
%             figure(5)
%             subplot(numPartitions,2,1)
%             figure(6)
%             subplot(numPartitions+1,2,1)
            [data] = decompAgglomProto3DProto3...
                                    (dataList{i},truncx,truncy,truncz,alphaList(i),voxelSizeList(i),YDFlag,dCalcFlag,dFile,dSaveFlag...
                                    ,numPartitions,dualLength,numIterations,relaxFlag,relaxIndex,precondFlag,pressureBasisFlag,upscalingFactor...
                                    ,agglomFlag,microFlag,agmgFlag,unixFlag,condDumpSerialiserFlag,diskDumpSerialiserFlag,preProcDiskDumpFlag,simulationDiskDumpFlag,runYDInSerial,runPartsInSerial,plotFlag);
%             solnList{i}=[data.timeVec,data.permVec,data.resVec];
%             Data=[timeVec,permVec,resVec];
            save([dataList{i},'Data',num2str(numPartitionsVec(j)),'PartsDL',num2str(dualLength),'Par.mat'],'data') 
        end
%     end
end

%%

% timeVec=data.timeVec;
% resVec=data.resVec;
% permVec=data.permVec;
% minPermVec=data.minPermVec;
% figure(3)
% semilogy(timeVec,resVec)
% drawnow
% figure(4)
% plot(timeVec,permVec);hold on;
% plot(timeVec,minPermVec);
% midPermVec=(minPermVec+permVec)/2;
% plot(timeVec,midPermVec);
% figure(5)
% plot(timeVec,permVec);hold on;
% plot(timeVec,minPermVec);
% midPermVec=(minPermVec+permVec)/2;
% plot(timeVec,midPermVec);




