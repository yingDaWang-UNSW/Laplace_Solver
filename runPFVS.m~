%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pore-Scale Finite Volume Solver + Agglomerated Dual Grid Domain Decomposition
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
cd /home/user/portableSoftwares/Dual_Grid_Laplace_Solver-master/
addpath(genpath(pwd))
mkdir transDump
%%
domainFile='geopack';
voxelSize=10*1e-6;
numPartitions=[2,2,2];
agglomFlag=1;

preProcDiskDumpFlag=1;
simulationDiskDumpFlag=0;
runPartsInSerial=0;
% clean out all disp lines
% 
% add velP output - careful to reduce memory usage - agglom doesnt yet
% support velp
[permeability,fields] = decompAgglom3D(domainFile,voxelSize,numPartitions,agglomFlag,preProcDiskDumpFlag,simulationDiskDumpFlag,runPartsInSerial);





