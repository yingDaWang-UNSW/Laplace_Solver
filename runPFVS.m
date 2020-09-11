%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pore-Scale Finite Volume Solver + Agglomerated Dual Grid Domain Decomposition
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
domainFile='geopack';
voxelSize=10*1e-6;
numPartitions=[2,2,2];
numIterations=50;
agglomFlag=1;

preProcDiskDumpFlag=0;
simulationDiskDumpFlag=0;
runPartsInSerial=0;

[data] = decompAgglom3D(domainFile,voxelSize,numPartitions,numIterations,agglomFlag,preProcDiskDumpFlag,simulationDiskDumpFlag,runPartsInSerial);





