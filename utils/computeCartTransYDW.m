function [T] = computeCartTransYDW(rock)
numCells=numel(rock.perm);
%get the trans  
cellNo = repmat(1:numCells,[6 1]);% rldecode(1 : numCells, diff(G.cells.facePos), 2) .';
cellNo = cellNo(:);

localFaces=[-0.5 0 0; 
            0.5 0 0;
            0 -0.5 0;
            0 0.5 0;
            0 0 -0.5;
            0 0 0.5];
C = repmat(localFaces,[numCells 1]);%G.cells.centroids;
% C = G.faces.centroids(G.cells.faces(:,1), :) - C(cellNo,:);

% cf  = G.cells.faces(:,1);
% sgn = 2*(cellNo == G.faces.neighbors(cf, 1)) - 1;
localNormals=[-1 0 0;
              1 0 0;
              0 -1 0;
              0 1 0;
              0 0 -1;
              0 0 1];
N = repmat(localNormals,[numCells 1]);%bsxfun(@times, sgn, G.faces.normals(cf, :));


K = rock.perm * [1, 0, 0, 1, 0, 1];

K = K(:, [1, 2, 3, 2, 4, 5, 3, 5, 6]);
i =      [1, 1, 1, 2, 2, 2, 3, 3, 3] ;
j =      [1, 2, 3, 1, 2, 3, 1, 2, 3] ;

T = zeros(size(cellNo));
for k = 1 : size(i, 2),
  T = T + (C(:, i(k)) .* K(cellNo, k) .* N(:, j(k)));
end
T = T ./ sum(C .* C, 2);

end