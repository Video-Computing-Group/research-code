function [bel, factorBels, cflag] = factor_gbp(localPot,factorPot,varargin)
% [bel, factorBels, cflag] = factor_gbp(localPot,factorPot,varagin)
% 
% Input arguments:
% 1. localPot: 1xN cell array, each cell is a column-vector, Vix1, defines the
%              local potentials for each node
% 2. factorPot: 1x<Num Factors> cell array, each cell {i} is a vector of
%               the potentials in factor i.
% 3. varargin: additional parameters needed for the algorithm:
%              'temperature', the temperature parameter
%                             the given local and pairwise potentials
%                             will be raised in power of 1/T
%                             default: 1.0
%              'regions',    a cell array with vectors of sorted nodes' indices,                
%                            where each such vector is a region (factor) in the graph
%                            default: each pair of neighbouring nodes is a
%                                     region (equals loopy)
%              'sum_or_max', the method for outgoing messages calculation
%                            0 - sum, 1 - max (default)
%              'alpha',      averaging parameter between new and old
%                            messages, where the new messages weight
%                            alpha
%                            default: 0.5
%              'trw',        0 - regular gbp, 1 - tree-rewighted gbp
%                            default: 0
%              'initMsg',    initial state for messages. should be
%                            1xNr cell array, where each cell {i} is
%                            a 1xnum-neighbours(i) cell array itself,
%                            such that {i}{n} contains a vector at length
%                            num-states(max(i,j)) (which means the
%                            number of states of the "son" region).
%                            this vector contains the messages from
%                            region i to region j.
%              'max_iter',   maximal number of iterations in inference
%                            default: 2000
%
%

SUM = 0;
MAX = 1;

% default values for specific-algorithm's arguments
gbp.sumOrMax = MAX;
gbp.alpha = 0.5;
gbp.regions = [];
gbp.all_regs = 0;
gbp.trw = 0;
gbp.counting_node = [];
gbp.initMsg = [];
% default value for temperature
temperature = 1.0;
% default value for max-iterations
maxIter = 2000;

adj = [];

% read given values for specific-algorithm's arguments
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
      case 'temperature'
          temperature = args{i+1};
      case 'sum_or_max'
          gbp.sumOrMax = args{i+1};
      case 'trw'
          gbp.trw = args{i+1};
      case 'initMsg'
          gbp.initMsg = args{i+1};
      case 'counting'
          gbp.counting_node = args{i+1};
      case 'alpha'
          gbp.alpha = args{i+1};
      case 'regions'
          gbp.regions = args{i+1};
      case 'adj'
          adj = args{i+1};
      case 'max_iter'
          maxIter = args{i+1};
      otherwise,
          error(['invalid argument name ' args{i}]);
  end
end

% if regions are not defined - pairs of neighbours are taken
% (+isolated nodes as seperate regions)
if (isempty(gbp.regions)==1)
    if (isempty(adj)==1)
        error('when the regions/factors are not given, the adjacencies matrix should be');
    else
        [I,J]=find(triu(adj));
        for i=1:length(I)
            gbp.regions{i} = [I(i), J(i)];
        end
        isolated = find(sum(adj)==0);
        for i=1:length(isolated)
            gbp.regions{end+1} = isolated(i);
        end
    end
end
if (isempty(gbp.counting_node)==1)
    gbp.counting_node = zeros(1,length(localPot));
    if (gbp.trw == 0)
        for i=1:length(gbp.regions)
            I = gbp.regions{i};
            gbp.counting_node(I) = gbp.counting_node(I) + 1;
        end
        gbp.counting_node = 1 - gbp.counting_node;
    end
end

[bel, factorBels, cflag] = c_factor_gbp(gbp.regions,localPot,factorPot,maxIter,...
    gbp.sumOrMax,gbp.alpha,gbp.trw,gbp.counting_node,gbp.initMsg,temperature);

