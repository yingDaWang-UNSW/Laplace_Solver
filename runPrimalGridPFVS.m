%% startup DGDD-PFVS
addpath(genpath(pwd))
%%

load geopack
domain=geopack;

voxelSize=1e-6;
alpha=1;
gradk=0;
FDGPA=1;
FDGPAYD=0;
micropore=0;
Pin=1;
Pout=0;
condFlag=0;
cond=[];
[K,vel,P,phi,phiEff]= solvePFVS(domain,voxelSize,alpha,gradk,FDGPA,FDGPAYD,micropore,Pin,Pout,condFlag,cond);