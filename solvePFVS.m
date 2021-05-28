% startup
% clc; clear;
% close all
function [K,vel,P,phi,phiEff]= solvePFVS(geoSmaller,voxelSize,alpha,gradk,FDGPA,FDGPAYD,micropore,Pin,Pout,condFlag,cond)
tic
%%
% startup;
% gravity reset off
% mrstModule add libgeometry incomp mpfa mimetic ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui ad-fi
% mrstModule add coarsegrid
% mrstModule add agmg
% % addpath(genpath('C:\Users\Ying\Dropbox\FVT2016\'))
% addpath(genpath(pwd))
%% import geosmall
% profile off
% profile -memory on;
disp('Loading Data')

% load BC.dat
% geoSmaller=reshape(BC,[80 30 30]);
% geoSmaller(:,:,1)=0;
% geoSmaller(:,:,end)=0;
% geoSmaller(:,1,:)=0;
% geoSmaller(:,end,:)=0;
% load('geocoal.mat');
% trunc=1:340;
% geoSmaller=geocoal(trunc,trunc,trunc);%(1:700,1:700,1:700);%geosmall;%(1:50,1:50,1:50);
% geoSmaller=permute(geoSmaller,[2 1 3]);
% clear bentheimer
% clear geopack
% voxelSize=10e-6;
% gradk=false;
% FDGPA=false;
% FDGPAYD=true;
% micropore=false;
hexGrid=false;
hexFile=false;
hexFileName='';
% geoSmaller=zeros(10,10,10);
% geoSmaller([1 10], :,:)=1;
% geoSmaller(:,[1 10],:)=1;
Nx=size(geoSmaller,1);
Ny=size(geoSmaller,2);
Nz=size(geoSmaller,3);
%% save to stl
%     for i=1:size(geoSmaller,3)
%         i
%         slice=squeeze(geoSmaller(:,:,i));
%         figure(1)
%         imagesc(slice);
%         drawnow
%         pause(0.1)
%     end
gridDimsX=cumsum(ones(size(geoSmaller,1),1));
gridDimsY=cumsum(ones(size(geoSmaller,2),1));
gridDimsZ=cumsum(ones(size(geoSmaller,3),1));


% [faces,vertices] = CONVERT_voxels_to_stl('BC.stl',true,~geoSmaller,gridDimsX,gridDimsY,gridDimsZ,'ascii');
phi=1-sum(sum(sum(geoSmaller)))/numel(geoSmaller);
%% generate laplace equation
if micropore
    disp('Generating full grid')
    G=cartGrid(size(geoSmaller));
    rock.perm=geoSmaller(:);
    rock.perm(rock.perm<0.5)=0;
    G=mcomputeGeometry(G);
    phiEff=phi;
%% remove disconnections, extract grainspace and porespace 
else 
    [geoSmaller,~,connFlag] = removeDisconnections2(geoSmaller);
    if ~connFlag
        K=0;
        sol=[];
        G=[];
        phiEff=0;
        return
    end
    phiEff=1-sum(sum(sum(geoSmaller)))/numel(geoSmaller);

    if ~hexGrid
        disp('Generating reduced grid from cartesian: Extracting pores and grains')
        rock.perm=geoSmaller(:);
        indexMap=[1:Nx*Ny*Nz]';
    %     %needed due to obnoxious memory redundancy associated with G
        G=generateReducedGridYDW(size(geoSmaller),indexMap(rock.perm==1|isnan(rock.perm)));
%         G = cartGrid(size(geoSmaller),size(geoSmaller));
%         active = indexMap(rock.perm==1|isnan(rock.perm));
%         G = extractSubgrid(G, active);
        clear indexMap
        G=mcomputeGeometry(G);
    else
        disp('Generating reduced grid from hexahedrons: Extracting pores and grains')
        nodes=[];
        elems=[];
        hexafile=[];
        if ~hexFile
            [nodes, elems] = processHex(geoSmaller);
        else
            hexafile=hexFileName;
        end
        G=hexahedralGridYDW(nodes,elems,hexafile);
    %     G=hexahedralGrid(nodes,elems);
        clear nodes elems
        G.cells.indexMap=[1:Nx*Ny*Nz]';
        rock.perm=geoSmaller(:);
        G.cells.indexMap(rock.perm==1|isnan(rock.perm))=[];
        G.cartDims=[Nx Ny Nz];
        
        [x, y, z]=meshgrid(1:Nx, 1:Ny, 1:Nz);
        
        G.cells.centroids=[x(:)-0.5 y(:)-0.5 z(:)-0.5];clear x y z
        G.cells.centroids=G.cells.centroids(G.cells.indexMap,:);
        G.cells.centroids=G.cells.centroids(:,[2,1,3]);
%         G=mcomputeGeometry(G);

        G.cells.volumes=ones(numel(G.cells.indexMap),1);
        %we can just use the top corner
        faceCoords1=G.nodes.coords(G.faces.nodes(1:4:end),:);
        faceCoords2=G.nodes.coords(G.faces.nodes(2:4:end),:);
        faceCoords3=G.nodes.coords(G.faces.nodes(3:4:end),:);
        faceCoords4=G.nodes.coords(G.faces.nodes(4:4:end),:);

        G.faces.centroids=(faceCoords1+faceCoords2+faceCoords3+faceCoords4)./4;
        clear faceCoords1 faceCoords2 faceCoords3 faceCoords4

    %     G=computeGeometryLiteYDW(G,'verbose',true);
%         plotCellData(G,ones(G.cells.num,1));
    end
    %remove cells with binary of 1
%     G=removeCells(G,G.cells.indexMap(rock.perm==1|isnan(rock.perm)));
%     rock.perm=rock.perm(rock.perm<1);
%     rock.perm=ceil(rock.perm)+1;
%     disp('Computing face centroids')
%    
    %% Find all cells that are connected to inflow/outflow (deprecated)
%     disp('Identifying flow path')
    
%     q   = ones(G.cells.num,1);
%     q   = processPartition(G, q);
%     qn  = unique(q);
%     disp(['Identified ',num2str(numel(qn)),' partitions within ',num2str(G.cells.num),' cells'])
%     [a,~] = histc(q,qn);
%     an=a>Nz;
%     qnn=qn(an);
%     disp(['Identified ',num2str(numel(qnn)),' large partitions within ',num2str(numel(qn)),' partitions'])
%     con = false(size(qnn));
%     rem = true(G.cells.num,1);
% %     temp=cell(size(qn));
%     disp('Analysing Partitions')
%     hwb = waitbar(0,['Analysing Partitions']);
%     for n=1:numel(qnn)
% %         n
%         i=qnn(n);
% %         if mod(i,1000)==1
% %             disp(['Analysing Partition ',num2str(i),' of ',num2str(numel(qn))])
% %         end
%       waitbar(n/numel(qnn),hwb);
%        ind = find(q==i);
%        if numel(ind)<Nz
% %            disp('Partition too small, moving on')
%            continue
%        end
%        cf = gridCellFaces(G,ind);
%        x  = G.faces.centroids(cf,3);
%        con(i) = (min(x)-eps <= 0) && (max(x)+eps >= Nz);
%        if con(i)
%           rem(ind) = false;
% %             temp{i}=ind;
%             disp(['Identified connected partition ',num2str(i),' of ',num2str(numel(qn))])
%        end
%     end
%     close(hwb)
% %     temp=cell2mat(temp);
%    disp(['Complete, ',num2str(sum(con)),' partitions identified'])
% %     figure
% %     plotCellData(G,rock.perm(~rem),~rem);
% %     plotGrid(G,rem,'FaceColor',.7*[1 1 1]);
% %     view(3);
%     clear q qn con cf x ind a hwb i n qnn trunc an
%     %%
%     if sum(rem)>0
%         %% Remove unconnected cells and initialize model
%         if sum(rem)==G.cells.num
%             disp('No connection at all, exiting')
%             return
%         end
%         G   = removeCells(G, rem);     
%         rock.perm    = rock.perm(~rem);
%         disp('Recomputing Geometry')
%         G=mcomputeGeometry(G);
%     else
%        disp('No cells removed, continuing') 
%     end
end

fluid     = initSingleFluid('mu',1,'rho',1);
resSol = initResSol(G, 0);
disp('Identifying wall faces')
[bFaces,bFaceCells] = boundaryFaces(G);
[~,inda,~]=unique(bFaceCells); clear bFaceCells
bFaces=bFaces(inda); clear inda
tol=1e-4;
% west = abs(G.faces.centroids(bFaces,3))<tol;
% east = abs(G.faces.centroids(bFaces,3)-Nz)<tol;
% north = abs(G.faces.centroids(bFaces,1))<tol;
% south = abs(G.faces.centroids(bFaces,1)-Nx)<tol;
% up = abs(G.faces.centroids(bFaces,2))<tol;
% down = abs(G.faces.centroids(bFaces,2)-Ny)<tol;
disp('Identifying boundary faces')
faceNumList=[1:G.faces.num]';
west = abs(G.faces.centroids(:,3))<tol;
east = abs(G.faces.centroids(:,3)-Nz)<tol;
north = abs(G.faces.centroids(:,1))<tol;
south = abs(G.faces.centroids(:,1)-Nx)<tol;
up = abs(G.faces.centroids(:,2))<tol;
down = abs(G.faces.centroids(:,2)-Ny)<tol;
% figure
% plotFaces(G,(east))
% hold on
% plotFaces(G,(west))
% hold on
% plotFaces(G,(north))
% hold on
% plotFaces(G,(south))
% hold on
% plotFaces(G,(up))
% hold on
% plotFaces(G,(down))


%%
if gradk
    disp('Calculating gradk conductivity')
    rock.perm=ones(numel(G.cells.indexMap),1);
    S=computeCartTransYDW(rock);%computeTrans(G,rock);
    bc=[];
    bc = addBC(bc, bFaces, 'pressure',0);
    clear north south up down
    src=[];
%     srcs=bwdistsc(geoSmaller);
%     srcs=srcs(:); 
%     srcs=srcs(srcs~=0);
%     src=addSource(src,[1:G.cells.num]',double(srcs./max(srcs(:))).*voxelSize^3);
    src=addSource(src,[1:G.cells.num]',1);

    % bc = addBC([], bFaces(east), 'pressure',0);
    r_psolve = @(state) incompTPFA(state, G, S, fluid,'bc',bc, 'src',src, 'LinSolve', @(A,b) agmg(A,b,[],[],[],1));
    resSol=r_psolve(resSol);
    dFluxCell=faceFlux2cellVelocity(G,resSol.flux);
    rock.perm=sum(abs(dFluxCell),2);
    alpha=rock.perm;
    alphaMap=zeros(G.cartDims);
    alphaMap(G.cells.indexMap)=alpha;
    alpha=alphaMap;
    gradkCoeff=(sum(abs(dFluxCell),2));
%     scaleFact=quantile(resSol.pressure,0.5)./quantile(((sum(abs(dFluxCell),2))),0.5);
%     rock.perm=resSol.pressure./voxelSize./scaleFact;

%     rock.perm=max(rock.perm)+min(rock.perm)-rock.perm;
%     figure(1)
%     plotCellData(G, rock.perm);view(3);axis equal
%     colorbar
%     kcart=zeros(Nz*Ny*Nx,1);
%     kcart(G.cells.indexMap)=rock.perm;
%     kcart=reshape(kcart,[Nx Ny Nz]);
%     for i=1:size(kcart,1)
%         slice=squeeze(kcart(:,:,i));
%         figure(1)
%         imagesc(slice);
%         drawnow
%         pause(1)
%     end
end
clear dFluxCell
%% scale pressure into conductivity FIX THIS SECTION, THERE ARE REDUNDANCIES
% alpha=2;
if ~condFlag
    if FDGPA
        disp('Calculating fdgpa conductivity')
        D=bwdist(geoSmaller);
        maxdist=max(D(:));
    % tic
        [Dmax] = FindWW2( D, maxdist,Nx,Ny,Nz, alpha,voxelSize );
        D=alpha./8.*(2.*Dmax.*D-D.^2);

    %     D=D./mean(D(:)).*gradkCoeff;
    %     D=D(:);
    %     figure(2)
    %     Dtest=D(D~=0);
    %     plotCellData(G, Dtest(:));view(3);axis equal
    %     colorbar

    % toc
    % tic
    elseif FDGPAYD
        disp('Calculating fdgpa conductivity')
        D=bwdist(geoSmaller);
        Dmax=calcDmax(D,24); 
    %     D=D./mean(D(:)).*gradkCoeff;
        D=alpha./8.*(2.*Dmax.*D-D.^2);
        
    else
        D=gradkCoeff;clear gradkCoeff
    end
else
    D=cond;
    clear cond
end

if ~micropore
    rock.perm=double(D(D~=0));
%     rock.perm    = rock.perm(~rem); clear rem
else
    D(D==0)=eps;
    rock.perm=double(D(D~=0));
end
% if gradk
% %     alpha=gradkCoeff./mean(rock.perm);
%     alpha=D./(rock.perm);
%     rock.perm =rock.perm.*alpha;
% end
clear D Dmax
% disp(['Conductivity calculated, alpha correction factor = ',num2str(mean(alpha))])
%% Dmax tester
%     rock.perm=double(Dmax(Dmax~=0));
%     rock.perm    = rock.perm(~rem); 
%     %%
%     Gcart=cartGrid(size(geoSmaller));
% 
% 
%     Dmaxlist=unique(Dmax);
%     Dmaxlist(Dmaxlist==0)=[];
%     for i=size(Dmaxlist,1):-1:1
%         tempvec=zeros(Nz*Ny*Nx,1);
%         tempvec(G.cells.indexMap)=rock.perm;
%         tempvec(tempvec~=Dmaxlist(i))=0;
%         figure(1)
%         plotCellData(Gcart,  tempvec,tempvec~=0);view(3);axis equal
%         drawnow;
%         hold on
%     end
% Dmax=Dmax./0.83;
% for i=1:size(Dmax,1)
%     slice=squeeze(Dmax(:,:,i));
%     figure(1)
%     imagesc(slice);
%     colorbar
%     
%     
%     slice=squeeze(Dmaxtrai(:,:,i));
%     figure(2)
%     imagesc(slice);
%     colorbar
% 
%     slice=squeeze(Dmaxtrai(:,:,i)-Dmax(:,:,i));
%     figure(3)
%     imagesc(slice);
%     colorbar
%     drawnow
% pause(1)
% end
% % toc
%%
clear geoSmaller
%% calculate AD
% N   = double(G.faces.neighbors);
% n    = size(N,1);
% neighbors = getNeighbourship(G, 'Topological', true);
% internal = all(neighbors~=0, 2);
% ic1  = neighbors(internal,1);
% ic2  = neighbors(internal,2);
% avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));

%% pass pressure to parm
disp('Calculating P and V')
% deltaFlux1=accumarray([ic1],[resSol.pressure],[G.cells.num 1]);
% deltaFlux2=accumarray([ic2],[resSol.pressure],[G.cells.num 1]);
% dk=deltaFlux1-deltaFlux2;
resSol = initResSol(G, 0);
bc=[];
% bc = addBC(bc, faceNumList(west), 'flux',Pin/numesl(faceNumList(west)));
bc = addBC(bc, faceNumList(west), 'pressure',Pin);

% bc = fluxside(bc,G,'Top',Pin);
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
% figure(2)
% plotCellData(G, resSol.pressure);view(3);axis equal
% colorbar
%% post process for avg perm
disp('Calculating K')
dFluxCell=faceFlux2cellVelocity(G,resSol.flux);
dFluxCell(isnan(dFluxCell)) = 0;

xvel2 = dFluxCell(:,1);%*(1/voxelSize);
yvel2 = dFluxCell(:,2);%*(1/voxelSize);
zvel2  = dFluxCell(:,3);%*(1/voxelSize);
clear dFluxCell
vMag = sqrt(xvel2.^2+yvel2.^2+zvel2.^2);

% dFluxCell=dFluxCell.*1e10
% vtkwrite('velocityfieldtemp.vtk','unstructured_grid',G.cells.centroids(:,1),G.cells.centroids(:,2),G.cells.centroids(:,3),...
%          'vectors','velocity',dFluxCell(:,1),dFluxCell(:,2),dFluxCell(:,3)); 
% writes a 3D unstructured grid that contains both vector and scalar values.
%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
%  positon, vector and scalar values.
%flow in z direction
Lz = Nz;
Lx = Nx;
Ly = Ny;

vMagcart=zeros(Lz*Ly*Lx,1);
vMagcart(G.cells.indexMap)=vMag;
vMagcart=reshape(vMagcart,[Lx Ly Lz]);

xVel2=zeros(Lz*Ly*Lx,1);
xVel2(G.cells.indexMap)=xvel2;
xVel2=reshape(xVel2,[Lx Ly Lz]);

yVel2=zeros(Lz*Ly*Lx,1);
yVel2(G.cells.indexMap)=yvel2;
yVel2=reshape(yVel2,[Lx Ly Lz]);

zVel2=zeros(Lz*Ly*Lx,1);
zVel2(G.cells.indexMap)=zvel2;
zVel2=reshape(zVel2,[Lx Ly Lz]);
% Area = (voxelSize)^2;
% q = zeros(Lz,1);
% tempuz = zeros(Lz,1);
% zVel2=zVel2.*Area;
% for k = 1:Lz
%     tempuz(k) = sum(sum(zVel2(:,:,k)));
% end
% q = tempuz(:);%.*Area;
% meanq = mean(q);
% K = -1.*(Lz/(Lx*Ly))*((meanq*1)/(Pin-Pout))/(voxelSize*1)*10^12;%Darcy;

vRMS2=(xVel2.^2+yVel2.^2+zVel2.^2).^0.5;
pressure=zeros(Lx,Ly, Lz);
pressure(G.cells.indexMap)=resSol.pressure;
Pin=pressure(:,:,1);
Pout=pressure(:,:,Lz);

Pin=mean(Pin(Pin>0));
Pout=mean(Pout(Pout>0));
gradP = (Pin-Pout)/Lz;

K = voxelSize^2*mean(vRMS2(:))/gradP*10^12;
Kzz = voxelSize^2*mean(zVel2(:))/gradP*10^12;
Kzx = voxelSize^2*mean(xVel2(:))/gradP*10^12;
Kzy = voxelSize^2*mean(yVel2(:))/gradP*10^12;

disp(['Permeability as calculated by PFVS: ',num2str(K), ' Darcies']);
sol=resSol;
K=[K, Kzx, Kzy, Kzz];
vel={xVel2,yVel2,zVel2};
P=pressure;
% [x,y,z]=meshgrid(1:Nx,1:Ny,1:Nz);
% % vecvtkydw(x,y,z,xvel2,yvel2,zvel2)
%  vtkwrite('velROI4down.vtk','structured_grid',x,y,z,'vectors','Velocity',xVel2,yVel2,zVel2);
%  writes
%  a structured 3D vector data into VTK file, with name specified by the string
%  filename. (u,v,w) are the vector components at the points (x,y,z). x,y,z
%  should be 3-D matrices like those generated by meshgrid, where
%  point(ijk) is specified by x(i,j,k), y(i,j,k) and z(i,j,k).
%  The matrices x,y,z,u,v,w must all be the same size and contain
%  corrresponding position and vector component. The string title specifies
%  the name of the vector field to be saved.

% figure(1)
% temp=faceFlux2cellVelocity(G,resSol.flux);
% temp=sqrt(sum(temp.^2,2));
% temp(isnan(temp))=0;
% G.nodes.coords=G.nodes.coords+repmat([0 0 localOrigin],[G.nodes.num 1]);
%             G=computeGeometryLiteYDW(G,'verbose',false);
%             plotFaceData(G,resSol.flux);view(3);axis equal;colorbar
%             plotCellData(G, temp);view(3);axis equal;colorbar
% plotCellData(G, resSol.pressure);view(3);axis equal;colorbar



% hold on
% drawnow
toc
% profile report
% profile off
% for i=1:size(zVel2,1)
%     slice=squeeze(zVel2(:,:,i));
%     figure(1)
%     imagesc(slice);
%     drawnow
%     pause(1)
% end
% plotCellData(G,vMag);view(3);axis equal

% plot(G.cells.centroids(:,3),zvel2,'o')
