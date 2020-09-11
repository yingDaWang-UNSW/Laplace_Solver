function [G] = generateReducedGridYDW(imageSize, cells)
%Remove cells from grid and renumber cells, faces and nodes.
G=cartGrid(imageSize);


   if isempty(cells),
      % Function is no-op if no cells are scheduled for removal.
      %
      cellmap = (1 : G.cells.num) .';
      facemap = (1 : G.faces.num) .';
      nodemap = (1 : G.nodes.num) .';

      return;
   end

   if islogical(cells),
      assert (numel(cells) == G.cells.num, ...
             ['The length of logical vector differs from ', ...
              'number of grid cells.']);
      cells = find(cells);
   end
  % New numbering of cells
  ind        = false(G.cells.num,1);
  ind(cells) = true;
% disp('Renumbering CellFaces')
  % remove and renumber cells in cellFaces
  if isfield(G.cells, 'numFaces'),
     warning('MRST:deprecated', ...
            ['The field numFaces will not be supported in ', ...
             'future releases of MRST.'])
     % 'numFaces' exists.  Preserve.
     numFaces                = G.cells.numFaces;
     G.cells.numFaces(cells) = [];
  else
     numFaces = diff(G.cells.facePos);
  end
  G.cells.faces(rldecode(ind, numFaces), :) = [];
% disp('Renumbering Neighbours')
  % Alter cell numbering in G.faces.neighbors
  n = G.faces.neighbors;
    cellmap    = mapExcluding(ind);
  G.faces.neighbors(n(:,1)>0,1) = cellmap(n(n(:,1)>0,1));
  G.faces.neighbors(n(:,2)>0,2) = cellmap(n(n(:,2)>0,2));
  clear n cellmap
% disp('Renumbering Cells')
  % Alter cells
  numFaces(cells)         = [];
  G.cells.num             = G.cells.num-numel(cells);
  G.cells.facePos         = cumsum([1; double(numFaces)]);
  clear numFaces
  if isfield(G.cells, 'indexMap'), G.cells.indexMap(cells) = []; end
%   disp('Altering Geometry')
  % Alter geometry
  if or(any(strcmp(G.type, 'computeGeometry')), any(strcmp(G.type, 'mcomputeGeometry')))
      G.cells.centroids(cells,:) = [];
      G.cells.volumes(cells,:)   = [];
      G.faces.areas(ind)         = [];
      G.faces.centroids(ind,:)   = [];
      G.faces.normals(ind,:)     = [];
  end
  clear cells
% disp('Renumbering Faces')
  %new numbering of faces.
  ind     = all(G.faces.neighbors(:,1:2)==0,2);

  % remove and renumber faces in faceNodes
  if isfield(G.faces, 'numNodes'),
     warning('MRST:deprecated', ...
            ['The field ''numNodes'' will be removed in a', ...
             'future release of MRST.'])
     % 'numNodes' exists.  Preserve.
     numNodes              = G.faces.numNodes;
     G.faces.numNodes(ind) = [];
  else
     numNodes = diff(G.faces.nodePos);
  end
  G.faces.nodes(rldecode(ind, numNodes)) = [];
  facemap = mapExcluding(ind);

  % remove and renumber faces in cellFaces
  G.cells.faces(:,1) = facemap(G.cells.faces(:,1));
clear facemap

  if any(G.cells.faces(:,1)==0),
      error('In removeCells: Too many faces removed!');
  end

  assert (isfield(G, 'type'), ...
          'Every grid must record its origin in the field ''type''.');


  numNodes(ind)            = [];
  G.faces.neighbors(ind,:) = [];
  G.faces.nodePos = cumsum([1; double(numNodes)]);
  clear numNodes
  if isfield(G.faces, 'tag'),
     G.faces.tag      (ind,:) = [];
  end
  G.faces.num              = G.faces.num - sum(ind);

%     disp('Removing Nodes')


  % Construct node map:
  ind = true(G.nodes.num, 1);
  ind(G.faces.nodes) = false;
  nodemap = mapExcluding(ind);

  % Remove nodes
  G.nodes.coords(ind,:) = [];
  G.nodes.num           = G.nodes.num - sum(ind);
  G.faces.nodes           = nodemap(G.faces.nodes);
clear nodemap
  if any(G.faces.nodes==0),
      error('In removeCells: Too many nodes removed!');
  end

  G.type = [G.type, { 'removeCells' }];

%   cellmap = find(cellmap > 0);
%   facemap = find(facemap > 0);
%   nodemap = find(nodemap > 0);
end
function m = mapExcluding(indices)
   n            = numel(indices);
   ind          = ones(n,1);
   ind(indices) = 0;
   m            = cumsum(ind);
   m(indices)   = 0;
end
