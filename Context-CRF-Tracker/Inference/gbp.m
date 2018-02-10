function [regMsg,bels,converged,belpairs] = gbp(regs,adj,regAdj,local,regPot,assignInd,bethe,extractSingle,extractPairs,varargin)
% [bels,converged,belpairs] = gbp(regs,adjMatrix,regAdjMatrix,Ls,regPot,assignInd,bethe,extractSingle,extractPairs,varargin)
% makes inference using gbp, after the preprocess was done
%
% Input arguments:
% 1. regs: 1xNr cell array, each cell {i} is a row vector with the indices
%          of nodes in the region i
% 2. adj: NnxNn matrix, defines neighbours in the graph
% 3. regAdj: NrxNr sparse matrix, defines the neighbours in the region graph
% 4. local: 1xNn cell array, each cell is a column-vector, Vix1, defines the
%           local potentials for the node i
% 5. regPot: 1xNr cell array, each cell {i} is a row vector of length Vr_i,
%            defines the local potentials for the region i
% 6. assignInd: 1xNr cell array, each cell {i} is a cell array of 1xneighbNum(i),
%               in which each cell {n} is:
%               1. empty if i>j - where j is the n'th neighbour of i
%               2. a column vector of length Vr_i with the indices of the
%                  corresponding assignment of j, if i<j
% 7. bethe: a row vector in length Nr, with the bethe of each region
%           (related to the double-count of the regions)
% 8. extractSingle: a row vector of length Nn, where each element (i) is the
%                   index of the region from which the belief for node i should
%                   be taken
% 9. extractPairs: 1xNn cell array, each cell {i} is a row vector of length
%                  neighbNum(i). each element (n) in this vector is:
%                  1. -1 if i>j - where j is the n'th neighbour of i
%                  2. the index of the region from which the pairwise beliefs
%                     of the pair <i,j> should be extracted, if i<j
% 10. varargin: additional parameters needed for gbp: 
%               'sum_or_max', the method for outgoing messages calculation
%                             0 - sum, 1 - max (default)
%               'alpha',      averaging parameter between new and old
%                             messages, where the new messages weight
%                             alpha
%                             default: 0.5
%               'max_iter',   maximal number of iterations in inference
%                             default: 2000
%               'initMsg',    initial messages for gbp. see description of
%                             output argument "regions' messages" for
%                             definition of structure of this argument.
%                             default: unique distribution
%      
% Note: Nn - number of nodes, Nr - number of regions in the region-graph
%       Vn - number of states for each node, Vr - number of states for
%       each region
%
% Output arguments:
% 1. regions' messages - 1xNr cell array, each cell {i} contains a
%                        1xneighbNum(i) cell array, where each cell {n} is
%                        a row vector 1xnumStates(max(i,j)) (which means
%                        the number of state of the "son" region). this
%                        vector contains the messages from region i to
%                        region j, to which gbp reached in the last
%                        iteration.
% 2. bel - 1xNn cell array, each cell {i} contains a column vector Vix1, where
%          bel{i}(xi) = EstimatedPr(i,xi)
% 3. convereged - the number of iterations took for the algorithm to converge,
%                 or -1 if the algorithm did not converge
% 4. belpairs - the resulted pairwise beliefs.
%


adjCell = adjMat2Cell(adj);
regsAdjCell = adjMat2Cell(regAdj);
N = length(adjCell);
Vn = zeros(1,N);
for i=1:N
    Vn(i) = length(local{i});
end
sumOrMax = 1;
alpha = 0.5;
maxIter = 2000;
initMsg = [];

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
      case 'sum_or_max'
          sumOrMax = args{i+1};
      case 'alpha'
          alpha = args{i+1};
      case 'max_iter'
          maxIter = args{i+1};
      case 'initMsg'
          initMsg = args{i+1};
      otherwise,
          error(['invalid argument name ' args{i}]);
  end
end

if (nargout > 3)
    [regMsg,bels,converged,belpairs_cform]=c_gbp(regs,adjCell,regsAdjCell,Vn,regPot,...
        assignInd,bethe,extractSingle,extractPairs,sumOrMax,alpha,maxIter,initMsg);
    belpairs = cform2SparseCell(belpairs_cform, adjCell);
else
    [regMsg,bels,converged]=c_gbp(regs,adjCell,regsAdjCell,Vn,regPot,assignInd,bethe,...
        extractSingle,extractPairs,sumOrMax,alpha,maxIter,initMsg);
end    
