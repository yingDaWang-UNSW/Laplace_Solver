function [resSol,dualPressures] = nonDirichletCorrectionSubroutine3D(G,S,NDFP,...
                                bc,localOrigin,localDual,plotFlag,iterNum,...
                                dualSideFlag,rmFlag,activeFaces,init,precond,agglomFlag)
    warning off
%     agglomFlag=false;
        %% for each subgrid, create mrst grid structure (cartesian full)
        if dualSideFlag || G.cells.num==0
            resSol=[];
            dualPressures=[];
            return
        end
                %% transform grid
        if ~agglomFlag %for now
            G.nodes.coords=G.nodes.coords+repmat([localOrigin],[G.nodes.num 1]);
        end
        faceIndexList=[1:G.faces.num]';
        %% perform sanity check on pressure locations
        for i=1:3
            for j=1:2
                if ~isempty(NDFP{i,j}) && numel(NDFP{i,j})>1
                    subFaceCentroids=G.faces.centroids(faceIndexList(G.faces.centroids(:,i)==bc(i,j)),:);
                    subNDFPCentroids=NDFP{i,j}(:,1:3)-localOrigin;
%                     assert(isequal(subFaceCentroids,subNDFPCentroids))
                    
                    % check the centroids match up and sort if they dont
%                     tempBCFaceInds=faceIndexList(bcFaces{i,j});
%                     tempBCFaceCentroids=G.faces.centroids(tempBCFaceInds,:);
                    if ~isequal(subFaceCentroids,subNDFPCentroids)
                        %check which rows mismatch
                        if i==1||i==2
                        NDFP{i,j}=sortrows(NDFP{i,j},3,'ascend');
                        else
                        NDFP{i,j}=sortrows(NDFP{i,j},2,'ascend');
                            
                        end
                        subNDFPCentroids=NDFP{i,j}(:,1:3)-localOrigin;
                        try
                            assert(isequal(subFaceCentroids,subNDFPCentroids))
                        catch
                            disp('WARNING!!! CENTROID VECTORS DONT MATCH UP.')
                            if size(subFaceCentroids,1)<size(subNDFPCentroids,1)
                                disp('Subsystem overspecified, attempting local rematch')
                                [~, inda, ~]=intersect(subNDFPCentroids, subFaceCentroids, 'stable', 'rows');
                                NDFP{i,j}=NDFP{i,j}(inda,:);
                            else
                                disp('Subsystem underspecified, attempting data interpolation.')
                                [diffCentroids,inda] = setdiff(subFaceCentroids, subNDFPCentroids, 'stable', 'rows');
                                inda=inda-(1:numel(inda))';
                                interpVals=repmat(mean(NDFP{i,j}(:,4)), [numel(inda),1]);%NDFP{i,j}()
                                diffCentroids=[diffCentroids, interpVals];
                                NDFP{i,j}=insertrows(NDFP{i,j},diffCentroids,inda);
                            end
                        end
                    end
%                     figure(103)
%                     scatter3(NDFP{i,j}(:,1),...
%                              NDFP{i,j}(:,2),...
%                              -1.*NDFP{i,j}(:,3),...
%                              ceil(NDFP{i,j}(:,4).*100),...
%                              NDFP{i,j}(:,4),'filled');colormap jet; axis equal,view(3)
%                          hold on
%                                            
%                     figure(104)
%                     scatter3(subFaceCentroids(:,1)+localOrigin(1),...
%                              subFaceCentroids(:,2)+localOrigin(2),...
%                              -1.*subFaceCentroids(:,3)-localOrigin(3),...
%                              ceil(NDFP{i,j}(:,4).*100),...
%                              NDFP{i,j}(:,4),'filled');colormap jet; axis equal,view(3)
%                          hold on
                end
                if  numel(NDFP{i,j})==1
                    NDFP{i,j}=[0 0 0 NDFP{i,j}];
                end
            end
        end
        
        %% use init if iterNum=1 and initial guess vector is procided instead of NDFP
        if iterNum==1
            initVec=init(:);
            if ~isempty(init)
                for i=1:3 %for each of the 6 faces
                    for j=1:2
                        if activeFaces(i,j) %if the face is active interpolate the bc pressure
                            subFaceCentroids=G.faces.centroids(faceIndexList(G.faces.centroids(:,i)==bc(i,j)),:);
                            NDFP{i,j}=subFaceCentroids;
                            % find rough locations of the bcs on the
                            % smaller grid
                            interpSubs=ceil(NDFP{i,j}./2);
                            interpSubs(interpSubs==0)=1;
                            interpInds=sub2ind(size(init),interpSubs(:,1),interpSubs(:,2),interpSubs(:,3));
                            NDFP{i,j}=[NDFP{i,j} initVec(interpInds)];
                        end
                        
                    end
                end
            end
        end

        %% define system and bcs
%         disp('Identifying boundary faces')
        BC=[];
        bcFaces=cell(3,2);
        src=[];
        for i=1:3 %for each of the 6 faces
            for j=1:2
                bcFaces{i,j}=G.faces.centroids(:,i)==bc(i,j); % find the centroids on the faces
                if activeFaces(i,j) %if the face is active, apply the dual bc presssure
%                     if numel(NDFP{i,j})>4
%                         % check the centroids match up
%                         tempBCFaceInds=faceIndexList(bcFaces{i,j});
%                         tempBCFaceCentroids=G.faces.centroids(tempBCFaceInds,:);
%                         if ~isequal(tempBCFaceCentroids,(NDFP{i,j}(:,1:3)-localOrigin))
%                             NDFP{i,j}=sortrows(NDFP{i,j},2,'ascend');
%                             assert(isequal(tempBCFaceCentroids,(NDFP{i,j}(:,1:3)-localOrigin)))
%                         end
%                     end
                    BC = addBC(BC, faceIndexList(bcFaces{i,j}), 'pressure',NDFP{i,j}(:,4));
                elseif i==3 %if inactive, but along the direction of flow, add the exterior bcs 
                    BC = addBC(BC, faceIndexList(bcFaces{i,j}), 'pressure',-1*(j-2));
%                     if j==1
% %                         BC = addBC(BC, faceIndexList(bcFaces{i,j}), 'pressure',1);
%                         src=addSource(src,1,1);
%                     elseif j==2
%                         src=addSource(src,1000,-1);
%                     end

                end
            end
        end
        if agglomFlag
%             BC = coarsenBC(G, BC); 
        end
%         zMinus = G.faces.centroids(:,3)==bc(3,1);
%         zPlus = G.faces.centroids(:,3)==bc(3,2);
%         yMinus = G.faces.centroids(:,2)==bc(2,1);
%         yPlus= G.faces.centroids(:,2)==bc(2,2);
%         xMinus = G.faces.centroids(:,1)==bc(1,1);
%         xPlus= G.faces.centroids(:,1)==bc(1,2);
%%
%         bc=[];
%         if activeFaces(3,1)
%             bc = addBC(bc, faceIndexList(zMinus), 'pressure',NDFP{3,1});
%             %check for pressure mismatch
%             tempCentroids=G.faces.centroids(faceIndexList(zMinus),:);
%         else
%             bc = addBC(bc, faceIndexList(zMinus), 'pressure',1);
%         end
%         if activeFaces(3,2)
%             bc = addBC(bc, faceIndexList(zPlus), 'pressure',NDFP{3,2});
%         else
%             bc = addBC(bc, faceIndexList(zPlus), 'pressure',0);
%         end
% %         if mod(iterNum,2)
%         if activeFaces(1,2)
%             bc = addBC(bc, faceIndexList(xPlus), 'pressure',NDFP{1,2});
%         end
%         if activeFaces(1,1)
%             bc = addBC(bc, faceIndexList(xMinus), 'pressure',NDFP{1,1});
%         end
%         if activeFaces(2,2)
%             bc = addBC(bc, faceIndexList(yPlus), 'pressure',NDFP{2,2});
%         end
%         if activeFaces(2,1)
%             bc = addBC(bc, faceIndexList(yMinus), 'pressure',NDFP{2,1});
%         end
%                             figure(105)
%                          
%                                      scatter3(G.faces.centroids(BC.face,1)+localOrigin(1),...
%                      G.faces.centroids(BC.face,2)+localOrigin(2),...
%                      -1.*G.faces.centroids(BC.face,3)-localOrigin(3),100.*ones(numel(BC.face),1),...
%                              BC.value,'filled');colormap jet; axis equal,view(3)
%                          hold on
%         end
        resSol = initResSol(G, 1);
        fluid     = initSingleFluid('mu',1,'rho',1);
%         S=computeTrans(G,rock);(A,b,restart,tol,maxit,verbose,x0,ijob)
        restart=10;
        tol=1e-6;
        maxit=100;
        verbose=0;
        if rmFlag
            r_psolve = @(state) incompTPFA(state, G, S, fluid,'bc',BC,'src',src);
        else
            r_psolve = @(state) incompTPFA(state, G, S, fluid,'bc',BC, 'src',src,'LinSolve', @(A,b) agmg(A,b,restart,tol,maxit,verbose,precond));
        end
        %% solve system
%         warning off
        resSol=r_psolve(resSol);
        warning on
%         disp(['Preconditioner Efficiency: ',num2str(1-norm(resSol.pressure-x0))])
        %% if grids are agglomerated, extract the entire pressure mapping 
        if agglomFlag
            dualPressures=[G.faces.centroids+localOrigin, resSol.facePressure];
            
        end
        %% for cartesian routine, this works. Extract mid pressures, only centroid and pressures, no overlap
        if ~dualSideFlag && ~agglomFlag
            dualPressures=struct();
            dualFaceInds=(G.faces.centroids(:,1)==localDual(1,1))|...
                         (G.faces.centroids(:,2)==localDual(2,1))|...
                         (G.faces.centroids(:,3)==localDual(3,1));

            dualFaceInds=faceIndexList(dualFaceInds);

            dualCentreLeftCells=G.faces.neighbors(dualFaceInds,1);
            dualCentreLeftFaces=dualFaceInds(dualCentreLeftCells~=0);
            
            dualCentreRightCells=G.faces.neighbors(dualFaceInds,2);
            dualCentreRightFaces=dualFaceInds(dualCentreRightCells~=0);

            dualPressures.left.pressure=resSol.facePressure(dualCentreLeftFaces);
            dualPressures.right.pressure=resSol.facePressure(dualCentreRightFaces);
            dualPressures.left.centroids=G.faces.centroids(dualCentreLeftFaces,:)+localOrigin;
            dualPressures.right.centroids=G.faces.centroids(dualCentreRightFaces,:)+localOrigin;
%             dualPressures.centroids=G.faces.centroids(dualFaceInds,:)+localOrigin;
%             dualPressures.pressure=resSol.facePressure(dualFaceInds);


%             figure(99)
% %             G.nodes.coords=G.nodes.coords+repmat([localOrigin],[G.nodes.num 1]);
%             plotFaces(G,dualFaceInds,resSol.facePressure(dualFaceInds)); colormap jet; axis equal,view(3)
%             hold on
% %             plotFaces(G,(dualCentreLeftFaces),'r','FaceAlpha',.3);
% %             plotFaces(G,(dualCentreRightFaces),'b','FaceAlpha',.3);
%             figure(100)
%             plotFaces(G,BC.face,BC.value); colormap jet; axis equal,view(3)
%             hold on
% % %                     end
%             figure(101)
%             plotFaces(G,dualCentreRightFaces,resSol.facePressure(dualCentreRightFaces)); colormap jet; axis equal,view(3)
%             plotFaces(G,dualCentreLeftFaces,resSol.facePressure(dualCentreLeftFaces)); colormap jet; axis equal,view(3)
%             hold on
%             
%             figure(102)
%             scatter3(dualPressures.left.centroids(:,1),...
%                      dualPressures.left.centroids(:,2),...
%                      -1.*dualPressures.left.centroids(:,3),...
%                      dualPressures.left.pressure.*100,...
%                      dualPressures.left.pressure,'filled');colormap jet; axis equal,view(3)
%             hold on
%             scatter3(dualPressures.right.centroids(:,1),...
%                      dualPressures.right.centroids(:,2),...
%                      -1.*dualPressures.right.centroids(:,3),...
%                      dualPressures.right.pressure.*100,...
%                      dualPressures.right.pressure,'filled');colormap jet; axis equal,view(3)
            % fig 99, 101, and 102 shoud match. they do, so not an issue with
            % mrst indexing. all match, so centroid and pressures are
            % correctly placed. somehow, the reinput bc from the last
            % setare scraambled though.
        end
        
        %% sanity check on face pressures out and input NDFP boundary conditions [see if facepressures are scrambled by ressol]
%         temp=resSol.facePressure(faceNumList(xMinus));
%         flag=(sum(NDFP{1,1}(:)-temp))==0
        %%
%         figure(2)
%         plotCellData(G, resSol.pressure);view(3);axis equal;colorbar
%         hold on
%         drawnow
        %% save face boundary system
%         pressure=resSol.pressure;
%         flux=resSol.flux
        
%         subPress{i}=resSol.pressure;
%         subVel{i}=faceFlux2cellVelocity(G,resSol.flux);
%         subBoundFaceFlux{i}=[resSol.flux(faceNumList(east)) resSol.flux(faceNumList(west))];
%         subBoundFacePres{i}=[resSol.facePressure(faceNumList(east)) resSol.facePressure(faceNumList(west))];
%         subDualFacePres{i}=[resSol.facePressure(faceNumList(dualCentre))];

    %% extract pressures
%     subBoundFaceFlux=cell2mat(subBoundFaceFlux);
%     subBoundFacePres=cell2mat(subBoundFacePres);
%     subDualFacePres=cell2mat(subDualFacePres);

%     subBoundFaceFlux=subBoundFaceFlux(:,2:end-1);
%     subBoundFacePres=subBoundFacePres(:,2:end-1);

%     subBoundFaceFlux(:,2:2:end)=subBoundFaceFlux(:,2:2:end).*-1;
%     flux=subBoundFaceFlux;%(sum(subBoundFaceFlux,2));
%     pressure=[subBoundFacePres(:,1) subDualFacePres subBoundFacePres(:,end)];
%     residual = abs(sum(diff(sum(flux))))/max(flux(:));
%     dFluxCell=cell2mat(subVel');
%     dFluxCell(isnan(dFluxCell)) = 0;
%     xvel2 = (1/voxelSize)*dFluxCell(:,1);
%     yvel2 = (1/voxelSize)*dFluxCell(:,2);
%     zvel2  = (1/voxelSize)*dFluxCell(:,3);
%     vMag = sqrt(xvel2.^2+yvel2.^2+zvel2.^2);
% 
% 
%     %flow in z direction
%     Lz = Nz;
%     Lx = Nx;
%     Ly = Ny;
% 
%     vMagcart=zeros(Lz*Ly*Lx,1);
%     vMagcart=vMag;
%     vMagcart=reshape(vMagcart,[Lx Ly Lz]);
% 
% 
%     xVel2=zeros(Lz*Ly*Lx,1);
%     xVel2=xvel2;
%     xVel2=reshape(xVel2,[Lx Ly Lz]);
% 
%     yVel2=zeros(Lz*Ly*Lx,1);
%     yVel2=yvel2;
%     yVel2=reshape(yVel2,[Lx Ly Lz]);
% 
%     zVel2=zeros(Lz*Ly*Lx,1);
%     zVel2=zvel2;
%     zVel2=reshape(zVel2,[Lx Ly Lz]);
%     Area = (voxelSize)^2*ones(Lz,1);
%     q = zeros(Lz,1);
%     tempuz = zeros(Lz,1);
%     for k = 1:Lz;
%         tempuz(k) = sum(sum(zVel2(:,:,k)));
%     end
%     q = tempuz(:).*Area;
%     meanq = mean(q);
%     K = (Lz/(Lx*Ly))*((meanq*1)/(1-0))/(voxelSize*1)*10^12;%Darcy;
%     disp(['Calculated perm: ',num2str(K),' calculated residual: ', num2str(residual), ' dual cycle iteration number: ', num2str(n)])

end