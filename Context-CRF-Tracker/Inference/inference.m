function [bel, converged, belpairs, msgs] = inference(adj,lambda,local,algorithm_name,varargin)
% [bel, convereged] = inference(adj,lambda,local,algorithm_name,varargin)
% makes inference in the given undirected graph with N nodes,
% and Vi possible values for each node i
% 
% Input arguments:
% 1. adj: NxN matrix, defines neighbours in the graph
% 2. lambda: there are 2 forms for lambda:
%            a. in general MRF algorithms (loopy, gbp, gibbs, mean-field) :
%               lambda should be a NxN cell array, defines the potentials
%               between the nodes. each cell {i,j} is a VixVj matrix,
%               where Vi is the number of states for node i, and Vj is for node j.
%            b. in PottsMRF alhorithms (monte-carlo algorithms which are planned
%               for Potts model, i.e. metropolis and the cluster algorithms
%               wolff and swendsen-wang) :
%               here lambda should be a NxN matrix, defines the strength of
%               the interactions between the nodes. each element (i,j) is a
%               scalar.
%               if Jij = lambda(i,j), and Psi{i,j} is the general MRF
%               matrix with the potentials between i and j, than:
%               Psi{i,j} = [exp(Jij), 1; 1, exp(Jij)];
% 3. local: 1xN cell array, each cell is a column-vector, Vix1, defines the
%           local potentials for each node.
% 4. algorithm_name: the name of algorithm to use for inference. possible values:
%                    'loopy', 'gbp', 'gibbs', 'wolff', 'swendsen-wang',
%                    'metropolis', 'mean-field'
% 5. varargin: additional parameters needed for the different algorithms
%              should be given in pairs, like:
%              'initial_state', [1 0 0 1 1]
%              I) for all algorithms:
%                 'temperature', the temperature parameter
%                                the given local and pairwise potentials
%                                will be raised in power of 1/T
%                                default: 1.0
%              II) for each different algorithm:
%                  a. for loopy:
%                     'sum_or_max', the method for outgoing messages calculation
%                                   0 - sum, 1 - max (default)
%                     'strategy',   the updating strategy
%                                   0 - sequential (default), 1 - parallel
%                     'trw',        running Tree-Reweighted-BP. the given
%                                   argument here can be either a matrix
%                                   where each element <i,j> in the matrix
%                                   is the rho-value for the edge between
%                                   node i and node j,
%                                   or a scalar, if all edges should have
%                                   the same rho-value.
%                                   default: regular BP (not TRBP)
%                     'log',        running loopy using log-space.
%                                   0 - regular BP (not using
%                                   log-space) is the default
%                                   1 - using log-space
%                                   if the given argument is 1, than the
%                                   temperature value is used (either given
%                                   by the user or set to 1)
%
%                                   IMPORTANT NOTE:
%                                   ---------------
%                                   when using log space, the given local
%                                   and pairwise potentials MUST BE THE
%                                   COST VALUES AND NOT THE ACTUAL
%                                   POTENTIALS. (where
%                                   cost_i(xi)=exp(-local_i(xi) /
%                                   temperature)
%
%                     'log_bels',   relevant when running loopy using log-space.
%                                   0 - return regular beliefs (not in
%                                   log-space) is the default
%                                   1 - return the beliefs in log-space
%                                   the connection between the 2 options is
%                                   bels = exp( - (log-bels) / temperature)
%
%                     'save_time',  it is possible to run loopy in a way
%                                   that wastes more space but saves
%                                   running time.
%                                   0 - no time saving (=saves space)
%                                     - this is the default
%                                   1 - use time saving (=use more space)
%
%                     'initMsg',    initial state for messages. should be
%                                   1xN cell array, where each cell {i} is
%                                   a 1xN cell array itself, such that
%                                   {i}{j} contains a vector at length
%                                   num-states(j), with the messages from
%                                   node i to node j
%                  b. for gbp:
%                     'regions',    a cell array with vectors of sorted nodes' indices,                
%                                   where each such vector is a region in the graph
%                                   default: each pair of neighbouring nodes is a
%                                            region (equals loopy)
%                                   IMPORTANT: when the next argument,
%                                   'all_regs', is 1, the 'regions'
%                                   argument must be a cell array of the
%                                   layers in the graph, size 1xnum_layers.
%                                   each layer is a cell array with
%                                   regions.
%                     'all_regs',   indicates wether the regions given are
%                                   only the big regions, or all regions
%                                   0 - only big regions (default)
%                                   1 - all regions
%                     'sum_or_max', the method for outgoing messages calculation
%                                   0 - sum, 1 - max (default)
%                     'alpha',      averaging parameter between new and old
%                                   messages, where the new messages weight
%                                   alpha
%                                   default: 0.5
%                     'trw',        0 - regular gbp, 1 - tree-rewighted gbp
%                                   default: 0
%                     'initMsg',    initial state for messages. should be
%                                   1xNr cell array, where each cell {i} is
%                                   a 1xnum-neighbours(i) cell array itself,
%                                   such that {i}{n} contains a vector at length
%                                   num-states(max(i,j)) (which means the
%                                   number of states of the "son" region).
%                                   this vector contains the messages from
%                                   region i to region j.
%                  c. for monte-carlo algorithms (gibbs, wolff, swendsen-wang):
%                     'initial_state',     row vector with initial state for sampling
%                                          default: choose randomly legal state for
%                                                   each node
%                     'burning_time',      num of transitions to wait before sampling
%                                          default: 100
%                     'sampling_interval', num of transitions to wait between samplings
%                                          default: 20
%                     'num_samples',       num of samples to take for estimation
%                                          default: 1000
%                  d. for all other algorithms (loopy, gbp, mean-field):
%                     'max_iter',   maximal number of iterations in inference
%                                   default: 2000
%
% Output arguments:
% bel - 1xN cell array, each cell {i} contains a column vector Vix1, where
%       bel{i}(xi) = EstimatedPr(i,xi)
% convereged - the number of iterations took for the algorithm to converge,
%              or -1 if the algorithm did not converge
% belpairs - relevant for loopy and gbp, the resulted pairwise beliefs.
% msgs - relevant for loopy and gbp, the messages reached by inference
%
% Examples:
% inference(adj,lambda,local,'gibbs','initial_state',startX,'sampling_interval',200)
% inference(adj,lambda,local,'loopy','sum_or_max',0)
%

LOOPY = 1;
GBP = 2;
GIBBS = 3;
WOLFF = 4;
SWENDSEN_WANG = 5;
METROPOLIS = 6;
MEAN_FIELD = 7;

SUM = 0;
MAX = 1;

SEQUENTIAL = 0;
PARALLEL = 1;

algorithm_names_array = {'loopy', 'gbp', 'gibbs', 'wolff', 'swendsen-wang', 'metropolis', 'mean-field'};
algo_num = length(algorithm_names_array);
algorithm = 0;
monte_carlo.chosen = 0;
for i=1:algo_num
    if (length(algorithm_name)==length(algorithm_names_array{i}))
        if (algorithm_name == algorithm_names_array{i})
            algorithm = i;
        end
    end
end
if (algorithm==0)
    error(['invalid algorithm ' algorithm_name]);
end

model = ~iscell(lambda); % 1 for potts model, 0 for general
if (~model && (algorithm==WOLFF || algorithm==SWENDSEN_WANG))
    error([algorithm_name ' works for potts model, please insert lambda as described in header']);
end

% default values for specific-algorithm's arguments
switch algorithm,
    case LOOPY
        loopy.sumOrMax = MAX;
        loopy.strategy = SEQUENTIAL;
        loopy.trw = 0;
        loopy.rho = 0;
        loopy.initMsg = [];
        loopy.saveTime = 0;
        loopy.counting_node = [];
case GBP
        gbp.sumOrMax = MAX;
        gbp.alpha = 0.5;
        gbp.regions = [];
        gbp.all_regs = 0;
        gbp.reg_pot = [];
        gbp.trw = 0;
	gbp.full = 0;
        gbp.counting_node = [];
        gbp.initMsg = [];
    case {GIBBS, WOLFF, SWENDSEN_WANG, METROPOLIS}
        monte_carlo.chosen = 1;
        monte_carlo.initial_state_given = 0;
        monte_carlo.burning_time = 100;
        monte_carlo.sampling_interval = 20;
        monte_carlo.num_samples = 1000;
end
% default value for temperature
temperature = 1.0;
% default value for max-iterations
maxIter = 2000;
% default value for threshold for convergence (log will be taken if
% log-space is asked)
th = 1e-8;
% log-space parameters
logsp = 0;
logbels = 0;

% read given values for specific-algorithm's arguments
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
      case 'temperature'
          temperature = args{i+1};
      case 'sum_or_max'
          if (algorithm~=LOOPY)
              if (algorithm==GBP)
                  gbp.sumOrMax = args{i+1};
              else
                  warning(['sum_or_max argument irrelevant for algorithm ' algorithm_name]);
              end
          else
              loopy.sumOrMax = args{i+1};
          end
      case 'trw'
          if (algorithm~=LOOPY)
              if (algorithm~=GBP)
                  warning(['trw argument irrelevant for algorithm ' algorithm_name]);
              else
                  gbp.trw = args{i+1};
              end
          else
              loopy.rho = rhoMat(args{i+1},adj);
              loopy.trw = 1;
          end
      case 'initMsg'
          if (algorithm~=LOOPY)
              if (algorithm~=GBP)
                  warning(['initMsg argument irrelevant for algorithm ' algorithm_name]);
              else
                  gbp.initMsg = args{i+1};
              end
          else
              loopy.initMsg = args{i+1};
          end
      case 'strategy'
          if (algorithm~=LOOPY)
              warning(['strategy argument irrelevant for algorithm ' algorithm_name]);
          else
              loopy.strategy = args{i+1};
          end
      case 'log'
          if (algorithm~=LOOPY && algorithm~=MEAN_FIELD && algorithm~=GBP)
              warning(['log argument irrelevant for algorithm ' algorithm_name]);
          else
              logsp = args{i+1};
          end
      case 'log_bels'
          if (algorithm~=LOOPY && algorithm~=MEAN_FIELD && algorithm~=GBP)
              warning(['log_bels argument irrelevant for algorithm ' algorithm_name]);
          else
              logbels = args{i+1};
          end
      case 'save_time'
          if (algorithm~=LOOPY)
              warning(['save_time argument irrelevant for algorithm ' algorithm_name]);
          else
              loopy.saveTime = args{i+1};
          end
      case 'trees_per_node'
          if (algorithm~=GBP)
              warning(['trees_per_node parameter irrelevant for algorithm ' algorithm_name]);
          else
              gbp.counting_node = sum(adj,1) - (args{i+1}(:))';
          end
      case 'full'
          if (algorithm~=GBP)
              warning(['full parameter irrelevant for algorithm ' algorithm_name]);
          else
              gbp.full = args{i+1};
          end
      case 'counting_node'
          if (algorithm~=GBP)
              if (algorithm~=LOOPY)
                  warning(['counting_node parameter irrelevant for algorithm ' algorithm_name]);
              else
                  loopy.counting_node = (args{i+1}(:))';
              end
          else
              gbp.counting_node = (args{i+1}(:))';
          end
      case 'alpha'
          if (algorithm~=GBP)
              warning(['alpha parameter irrelevant for algorithm ' algorithm_name]);
          else
              gbp.alpha = args{i+1};
          end
      case 'regions'
          if (algorithm~=GBP)
              warning(['regions argument irrelevant for algorithm ' algorithm_name]);
          else
              gbp.regions = args{i+1};
          end
      case 'all_regs'
          if (algorithm~=GBP)
              warning(['all_regs argument irrelevant for algorithm ' algorithm_name]);
          else
              gbp.all_regs = args{i+1};
          end
      case 'regions_potentials'
          if (algorithm~=GBP)
              warning(['regions_potentials argument irrelevant for algorithm ' algorithm_name]);
          else
              gbp.reg_pot = args{i+1};
          end
      case 'initial_state'
          if (monte_carlo.chosen == 0)
              warning(['initial_state argument irrelevant for algorithm ' algorithm_name]);
          else
              monte_carlo.initial_state = args{i+1};
              monte_carlo.initial_state_given = 1;
          end
      case 'burning_time'
          if (monte_carlo.chosen == 0)
              warning(['burning_time argument irrelevant for algorithm ' algorithm_name]);
          else
              monte_carlo.burning_time = args{i+1};
          end
      case 'sampling_interval'
          if (monte_carlo.chosen == 0)
              warning(['sampling_interval argument irrelevant for algorithm ' algorithm_name]);
          else
              monte_carlo.sampling_interval = args{i+1};
          end
      case 'num_samples'
          if (monte_carlo.chosen == 0)
              warning(['num_samples argument irrelevant for algorithm ' algorithm_name]);
          else
              monte_carlo.num_samples = args{i+1};
          end
      case 'max_iter'
          if (monte_carlo.chosen == 1)
              warning(['max_iter argument irrelevant for algorithm ' algorithm_name]);
          else
              maxIter = args{i+1};
          end
      case 'threshold'
          if (monte_carlo.chosen == 1)
              warning(['max_iter argument irrelevant for algorithm ' algorithm_name]);
          else
              th = args{i+1};
          end
      otherwise,
          error(['invalid argument name ' args{i}]);
  end
end

if (logsp)
    th = log(th);
end

adjCell = adjMat2Cell(adj);
psiCell = [];
if ~isempty(lambda)
	psiCell = psiSqueezeCell(lambda,adjCell,model);
end
clear lambda;
if isa(local,'sparse_cell')
    local = full(local);
end

belpairs = [];
switch algorithm,
    case LOOPY
        if (logsp == 0 && loopy.trw ~= 0)
            epsilon = 10^(-64);
            for i=1:length(psiCell)
                if iscell(psiCell{i})
                    for n=1:length(psiCell{i})
                        psiCell{i}{n} = psiCell{i}{n} / max(max(psiCell{i}{n}));
                        psiCell{i}{n} = psiCell{i}{n} .^ (1/loopy.rho{i}(n));
                        psiCell{i}{n}(find(psiCell{i}{n}<epsilon)) = epsilon;
                    end
                else
                    psiCell{i} = psiCell{i} ./ loopy.rho{i};
                end
            end
        end
        clear adj;
        if (nargout > 3)
            [bel, converged, belpairs_cform, msgs] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,...
                maxIter,th,logsp,logbels,loopy.sumOrMax,loopy.strategy,loopy.trw,loopy.rho,loopy.initMsg,loopy.saveTime,loopy.counting_node);
            belpairs = cform2SparseCell(belpairs_cform, adjCell);
        else
            if (nargout > 2)
                [bel, converged, belpairs_cform] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,...
                    maxIter,th,logsp,logbels,loopy.sumOrMax,loopy.strategy,loopy.trw,loopy.rho,loopy.initMsg,loopy.saveTime,loopy.counting_node);
                belpairs = cform2SparseCell(belpairs_cform, adjCell);
            else
                [bel, converged] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,maxIter,th,...
                    logsp,logbels,loopy.sumOrMax,loopy.strategy,loopy.trw,loopy.rho,loopy.initMsg,loopy.saveTime,loopy.counting_node);
            end
        end
    case GBP
        if (isempty(gbp.counting_node)==1)
            if gbp.trw>0
                gbp.counting_node = zeros(1,size(adj,1));
            end
        else
            if (isempty(gbp.regions)==0)
                warning('since regions are defined, the trees_per_node parameter will be ignored');
                gbp.counting_node = [];
            end
        end
        % if regions are not defined - pairs of neighbours are taken
        % (+isolated nodes as seperate regions)
        if (isempty(gbp.regions)==1)
            [I,J]=find(triu(adj));
            for i=1:length(I)
                gbp.regions{i} = [I(i), J(i)];
            end
            isolated = find(sum(adj)==0);
            for i=1:length(isolated)
                gbp.regions{end+1} = isolated(i);
            end
        end

        gbp.full = gbp.full || ~gbp.trw;

        clear adj;
        if (nargout > 3)
            [bel, converged, belpairs_cform,msgs] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,maxIter,th, ...
                logsp,logbels,gbp.regions,gbp.all_regs,gbp.sumOrMax,gbp.alpha,gbp.trw,gbp.counting_node,gbp.full,...
		gbp.initMsg,gbp.reg_pot,0);

            belpairs = cform2SparseCell(belpairs_cform, adjCell);
        else
            if (nargout > 2)
                [bel, converged, belpairs_cform] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,maxIter,th, ...
                	logsp,logbels,gbp.regions,gbp.all_regs,gbp.sumOrMax,gbp.alpha,gbp.trw,gbp.counting_node,gbp.full,...
			gbp.initMsg,gbp.reg_pot,0);

                belpairs = cform2SparseCell(belpairs_cform, adjCell);
            else
                [bel, converged] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,maxIter,th, ...
                	logsp,logbels,gbp.regions,gbp.all_regs,gbp.sumOrMax,gbp.alpha,gbp.trw,gbp.counting_node,gbp.full,...
			gbp.initMsg,gbp.reg_pot,0);
            end
        end
    case MEAN_FIELD
        clear adj;
        [bel, converged] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model,maxIter,th,logsp,logbels);
    case {GIBBS, WOLFF, SWENDSEN_WANG, METROPOLIS}
        if (monte_carlo.initial_state_given == 0)
            num_nodes = length(local);
            monte_carlo.initial_state = zeros(1,num_nodes);
            for i=1:num_nodes
                vi = length(local{i});
                q = ones(1,vi) / vi;
                monte_carlo.initial_state(i) = chooseInteger(q);
            end
        end
        clear adj;
        [bel, converged] = c_inference(adjCell,psiCell,local,algorithm-1,temperature,model, ...
            monte_carlo.initial_state,monte_carlo.burning_time, ...
            monte_carlo.sampling_interval,monte_carlo.num_samples);
    otherwise,
        error(['invalid algorithm ' algorithm_name]);
end
