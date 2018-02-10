function [allRegs,regAdj,regPot,assignInd,bethe,extractSingle,extractPairs,regsBP] = gbp_preprocess(adj,lambda,local,varargin)
% [regs,regAdj,regPot,assignInd,bethe,extractSingle,extractPairs] = gbp_preprocess(adj,lambda,local,varargin)
% runs in c++ the preprocessor for gbp
%
% Input arguments:
% 1. adj: NnxNn matrix, defines neighbours in the graph
% 2. lambda: NnxNn cell array, defines the potentials between the nodes.
%            each cell {i,j} is a VixVj matrix, where Vn_i is the number of
%            states for node i, and Vj is for node j.
% 3. local: 1xNn cell array, each cell is a column-vector, Vix1, defines the
%           local potentials for the node i
% 4. regs: a cell array with vectors of sorted nodes' indices,
%          where each such vector is a region in the graph
%          default: each pair of neighbouring nodes is a region (equals loopy)
%          IMPORTANT: when the argument 'all_regs' is 1, this argument
%                     must be a cell array of the layers in the graph,
%                     size 1xnum_layers.
%                     each layer is a cell array with regions.
% 5. varagin: additional parameters
%             'all_regs': indicates wether the regions given are only the
%                         big regions, or all regions
%                         0 - only big regions (default)
%                         1 - all regions
%             'trw': indicates wether Tree-Reweighted should be used.
%              
%
% Output arguments:
% 1. allRegs: 1xNr cell array, each cell {i} is a row vector with the nodes in
%             the region i (in the full graph, including the big regions
%             given)
% 2. regAdj: NrxNr sparse matrix, defines the neighbours in the region graph
% 3. regPot: 1xNr cell array, each cell {i} is a row vector of length Vr_i,
%            defines the local potentials for the region i
% 4. assignInd: 1xNr cell array, each cell {i} is a cell array of 1xneighbNum(i),
%               in which each cell {n} is:
%               1. empty if i>j - where j is the n'th neighbour of i
%               2. a column vector of length Vr_i with the indices of the
%                  corresponding assignment of j, if i<j
% 5. bethe: a row vector in length Nr, with the bethe of each region
%           (related to the double-count of the regions)
% 6. extractSingle: a row vector of length Nn, where each element (i) is the
%                   index of the region from which the belief for node i should
%                   be taken
% 7. extractPairs: 1xNn cell array, each cell {i} is a row vector of length
%                  neighbNum(i). each element (n) in this vector is:
%                  1. -1 if i>j - where j is the n'th neighbour of i
%                  2. the index of the region from which the pairwise beliefs
%                     of the pair <i,j> should be extracted, if i<j
%
% Note: Nn - number of nodes, Nr - number of regions in the region-graph
%       Vn - number of states for each node, Vr - number of states for
%       each region
%
% Important: in "regAdj" the c-indices [0...Nr-1] are converted to
%            matlab-indices [1...Nr], but not in "assignInd",
%            "extractSingle" and "extractPairs".
%            all those 4 arguments should be passed to gbp as is, and should not
%            be changed by the user.
%            the user may change ONLY "regPot" and "bethe".
%            "allRegs" is returned as information for the user, and should not
%            be passed to gbp later. only the big regions, given here as
%            "regs" input argument, ahould be given to gbp later.
%

all_regs = 0;
trw = 0;
gfull = 0;
counting_node = [];
regs = [];

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
      case 'all_regs'
          all_regs = args{i+1};
      case 'trw'
          trw = args{i+1};
      case 'full'
          gfull = args{i+1};
      case 'trees_per_node'
          counting_node = sum(adj,1) - (args{i+1}(:))';
      case 'regions'
          regs = args{i+1};
  end
end

adjCell = adjMat2Cell(adj);
psiCell = psiSqueezeCell(lambda,adjCell);

if (isempty(counting_node)==1)
    if trw>0
        counting_node = zeros(1,size(adj,1));
    end
else
    if (isempty(regs)==0)
        warning('since regions are defined, the trees_per_node parameter will be ignored');
        counting_node = [];
    end
end

if (isempty(regs)==1) || (nargout>7)
    regsBP = {};
    [I,J]=find(triu(adj));
    for i=1:length(I)
        regsBP{i} = [I(i), J(i)];
    end
    isolated = find(sum(adj)==0);
    for i=1:length(isolated)
        regsBP{end+1} = isolated(i);
    end
end
if (isempty(regs)==1)
    regs = regsBP;
end

gfull = gfull || ~trw;

[allRegs,regAdjCell,regPot,assignInd,bethe,extractSingle,extractPairs]=c_gbp_preprocess(adjCell,psiCell,local,regs,all_regs,trw,counting_node,gfull);
regAdj = adjCell2Mat(regAdjCell);
