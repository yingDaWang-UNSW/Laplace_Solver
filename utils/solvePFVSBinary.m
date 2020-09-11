function [K,sol,G,phi,phiEff,Pin,Pout]= solvePFVSBinary(geoSmaller,voxelSize,alpha,Pin,Pout)
tic
disp('Loading Data')
Nx=size(geoSmaller,1);
Ny=size(geoSmaller,2);
Nz=size(geoSmaller,3);
phi=1-sum(sum(sum(geoSmaller)))/numel(geoSmaller);
%% generate laplace equation
[geoSmaller,~,connFlag] = removeDisconnections(geoSmaller);
if ~connFlag
    K=0;
    sol=[];
    G=[];
    phiEff=0;
    return
end
phiEff=1-sum(sum(sum(geoSmaller)))/numel(geoSmaller);
disp('Generating reduced grid from cartesian: Extracting pores and grains')
rock.perm=geoSmaller(:);
indexMap=[1:Nx*Ny*Nz]';
G=generateReducedGridYDW(size(geoSmaller),indexMap(rock.perm==1|isnan(rock.perm)));
clear indexMap
[x, y, z]=meshgrid(1:Nx, 1:Ny, 1:Nz);
G.cells.centroids=[x(:)-0.5 y(:)-0.5 z(:)-0.5];clear x y z
G.cells.centroids=G.cells.centroids(G.cells.indexMap,:);
G.cells.centroids=G.cells.centroids(:,[2,1,3]);
G.cells.volumes=ones(numel(G.cells.indexMap),1);
faceCoords1=G.nodes.coords(G.faces.nodes(1:4:end),:);
faceCoords2=G.nodes.coords(G.faces.nodes(2:4:end),:);
faceCoords3=G.nodes.coords(G.faces.nodes(3:4:end),:);
faceCoords4=G.nodes.coords(G.faces.nodes(4:4:end),:);
G.faces.centroids=(faceCoords1+faceCoords2+faceCoords3+faceCoords4)./4;
clear faceCoords1 faceCoords2 faceCoords3 faceCoords4
fluid = initSingleFluid('mu',1,'rho',1);
disp('Identifying wall faces')
tol=1e-4;

disp('Identifying boundary faces')
faceNumList=[1:G.faces.num]';
west = abs(G.faces.centroids(:,3))<tol;
east = abs(G.faces.centroids(:,3)-Nz)<tol;

disp('Calculating fdgpa conductivity')
D=bwdist(geoSmaller);
maxdist=max(D(:));
[Dmax] = FindWW2( D, maxdist,Nx,Ny,Nz,alpha,voxelSize );
D=alpha./8.*(2.*Dmax.*D-D.^2);
rock.perm=double(D(D~=0));
clear D Dmax
clear geoSmaller
%% pass pressure to parm
disp('Calculating P and V')
resSol = initResSol(G, 0);
bc=[];
bc = addBC(bc, faceNumList(west), 'pressure',Pin);
bc = addBC(bc, faceNumList(east), 'pressure',Pout);
if isempty(bc)
    K=0;
    sol=struct();
    return
end
clear west east
S=computeCartTransYDW(rock);%computeTrans(G,rock); clear rock
r_psolve = @(state) incompTPFA(state, G, S, fluid,'bc',bc, 'LinSolve', @(A,b)  agmg(A,b,[],[],[],1));
resSol=r_psolve(resSol);
%% post process for avg perm
disp('Calculating K')
dFluxCell=faceFlux2cellVelocity(G,resSol.flux);
dFluxCell(isnan(dFluxCell)) = 0;

xvel2 = dFluxCell(:,1);%*(1/voxelSize);
yvel2 = dFluxCell(:,2);%*(1/voxelSize);
zvel2  = dFluxCell(:,3);%*(1/voxelSize);
clear dFluxCell
vMag = sqrt(xvel2.^2+yvel2.^2+zvel2.^2);

vMagcart=zeros(Nz*Ny*Nx,1);
vMagcart(G.cells.indexMap)=vMag;
vMagcart=reshape(vMagcart,[Nx Ny Nz]);

xVel2=zeros(Nz*Ny*Nx,1);
xVel2(G.cells.indexMap)=xvel2;
xVel2=reshape(xVel2,[Nx Ny Nz]);

yVel2=zeros(Nz*Ny*Nx,1);
yVel2(G.cells.indexMap)=yvel2;
yVel2=reshape(yVel2,[Nx Ny Nz]);

zVel2=zeros(Nz*Ny*Nx,1);
zVel2(G.cells.indexMap)=zvel2;
zVel2=reshape(zVel2,[Nx Ny Nz]);

vRMS2=(xVel2.^2+yVel2.^2+zVel2.^2).^0.5;
pressure=zeros(Nx,Ny, Nz);
pressure(G.cells.indexMap)=resSol.pressure;
Pin=pressure(:,:,1);
Pout=pressure(:,:,Nz);

Pin=mean(Pin(Pin>0));
Pout=mean(Pout(Pout>0));
gradP = (Pin-Pout)/Nz;

K = voxelSize^2*mean(vRMS2(:))/gradP*10^12;
Kzz = voxelSize^2*mean(zVel2(:))/gradP*10^12;
Kzx = voxelSize^2*mean(xVel2(:))/gradP*10^12;
Kzy = voxelSize^2*mean(yVel2(:))/gradP*10^12;

disp(['Permeability as calculated by PFVS: ',num2str(K), ' Darcies']);
sol=resSol;
K=[K, Kzx, Kzy, Kzz];

toc

