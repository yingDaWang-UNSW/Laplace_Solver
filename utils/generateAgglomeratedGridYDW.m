function [CG, rock] = generateAgglomeratedGridYDW(bwD,D,Dmax,origin, R)
%% defining where to upscale
    % regerenate the binary geometry
    sample = D==0;
    Ax = size(sample,1);
    Ay = size(sample,2);
    Az = size(sample,3);
%     R = 10e-6 ;
    rho = 1;
    mu = 1;
    Shape = 1;
    %%
    level = 1;
    downfac = 1/(2^level);
    xDown = Ax*downfac;
    yDown = Ay*downfac;
    zDown = Az*downfac;
    DD = Dmax;
    geo = sample;
    Disdown=zeros(xDown,yDown,zDown);
    DmaxDown = zeros(xDown,yDown,zDown);
    sizelim = 1/downfac;
    for i=1:xDown
        for j=1:yDown
            for k=1:zDown
                x= sizelim*i-(sizelim-1):sizelim*i;
                y= sizelim*j-(sizelim-1):sizelim*j;
                z= sizelim*k-(sizelim-1):sizelim*k;
               locD=bwD(x,y,z);
               locDmax = DD(x,y,z);
               Fine = geo(x,y,z);
                if sum(Fine(:)) > 0
                    condArithLoc= 0;
                    avgDmax = 0;
                else
                    avgDmax = mode(locDmax(:));
                    condArithLoc=  min(locD(:)) + 0.5;
                end
                % distribute to DDown
                DmaxDown(i,j,k) = avgDmax;
                Disdown(i,j,k)= condArithLoc;
            end
        end
    end
    %%
    
    sample1 = Disdown;
    Ax1 = size(sample1,1);
    Ay1 = size(sample1,2);
    Az1 = size(sample1,3);
    DDown2=zeros(Ax1, Ay1, Az1);
    for i = 1:Ax1
        for j = 1:Ay1
            for k = 1:Az1
                if Disdown(i,j,k) ~= 0
                    DDown2(i,j,k) = (Shape) * R^3 * (rho/(8*mu)) * (2 * DmaxDown(i,j,k) * (Disdown(i,j,k)) - ((Disdown(i,j,k)))^2 ) ;
                end
            end
        end
    end
    %%
    disp('Generating Agglomerated Grid')
    G = cartGrid([Ax, Ay, Az]);
    active = find(sample == 0);
    G = extractSubgrid(G, active);
    clear active
    G = mcomputeGeometry(G);
    %%
%     % Retrieve Cartesian indices of each active cell.
    [i, j, k] = gridLogicalIndices(G);
%     % Identify Cartesian indices of coarse blocks.  Assumes 2-by-2-by-2 blocks.
%     blk = cellfun(@(ix) reshape(floor(((1 : max(ix(:))) - 1) ./ 2) + 1, [], 1), { i, j, k }, 'UniformOutput', false);
%     % Generate initial partition vector based on 2-by-2-by-2 groups of cells
%     dims = cellfun(@(ix) ix(end), blk);  % ix(end) == max(ix)
%     p = sub2ind(dims, blk{1}(i(:)), blk{2}(j(:)), blk{3}(k(:)));
%     % Identify those blocks that must be refined
%     refine = UpMap(p) == 1;  % Linear indexing, relies on ordering details of SUB2IND
%     % Assign new block IDs for refined blocks (one block per cell in refinement region)
%     p(refine) = max(p) + (1 : sum(refine)).';
%     % Prune empty blocks.
%     p = compressPartition(p);
%     CG = generateCoarseGrid(G, p);
%     CG = coarsenGeometry(CG);
%     CG.faces.nodePos = CG.faces.connPos;
%     CG.faces.nodes = CG.faces.fconn;
%     CG.parent.nodes.coords = CG.parent.nodes.coords+origin;
%     
%     CG.faces.areas(CG.faces.areas>1)=1;
%     CG.cells.volumes(CG.cells.volumes>1)=1;
    p = 1:size(i);
    p = p(:);
    level = 1;

    downfac = 1/(2^level);
    xDown = Ax*downfac;
    yDown = Ay*downfac;
    zDown = Az*downfac;
    disp('Generating Partition array')
    mergemap = zeros(Ax,Ay,Az);
    sizelim = 2;
    count = max(p)+5;
    for i=2:xDown-1
        for j=2:yDown-1
            for k = 2:zDown-1
                x= sizelim*i-(sizelim-1):sizelim*i;
                y= sizelim*j-(sizelim-1):sizelim*j;
                z= sizelim*k-(sizelim-1):sizelim*k;
                Fine = sample(x,y,z);
                if sum(Fine(:)) > 0 || ismember(i,[xDown/2, xDown/2+1])|| ismember(j,[yDown/2, yDown/2+1])|| ismember(k,[zDown/2, zDown/2+1])
                    mergemap(x,y,z) = 0;             
                else
                    count = count + 1;
                    mergemap(x,y,z)  = count;    
                end
            end
        end
    end

    mergemap = mergemap + (~sample);
    mergemap = mergemap(mergemap>0);
    mergemap = mergemap(:);
    mergemap(mergemap==1) =0;
    temp = mergemap(:) > 0;
    p(temp) = mergemap(temp);
    p = compressPartition(p);
    CG = generateCoarseGrid(G, p);
    CG = coarsenGeometry(CG);
%     figure(222);plotGrid(CG, 150000:200000);

    CG.cells.volumes(CG.cells.volumes>1) =1;
    CG.faces.areas(CG.faces.areas>1) =1;
    % close all
%     figure(115); plotGrid(CG);%plotCellData(CG.parent,p);
%     view(3); camproj perspective; axis tight;axis equal
%     faceInds=1:CG.faces.num;
%     figure(111); 
%     plotFaces(CG, faceInds(CG.faces.centroids(:,1)==Ax/2));view(3); camproj perspective; axis tight;hold on
%     plotFaces(CG, faceInds(CG.faces.centroids(:,2)==Ay/2));view(3); camproj perspective; axis tight;hold on
%     plotFaces(CG, faceInds(CG.faces.centroids(:,3)==Az/2));view(3); camproj perspective; axis tight;hold on
%     axis equal;drawnow
    
    CG.parent=[];

    %% Assigning conductivity
    disp('Assigning Agglomerated conductivity values')
    rock.perm=ones(CG.cells.num,1);% = makeRock(CG, 1, 1);
    AA = size(rock.perm,1);


    for xx = 1:AA

      temp =  CG.cells.centroids(xx,:);

      xi = round(temp(1),1);
      yi = round(temp(2),1);
      zi = round(temp(3),1);

      if mod(xi,1) ==0.5 || mod(yi,1) ==0.5 || mod(zi,1)==0.5

          i = round(xi+0.5);
          j = round(yi+0.5);
          k = round(zi+0.5);

          rock.perm(xx) = D(i,j,k);

      else

          i = round(xi/2 + 0.5);
          j = round(yi/2 + 0.5);
          k = round(zi/2 + 0.5);
          rock.perm(xx) = DDown2(i,j,k);
      end
    end

    clear D DDown
end