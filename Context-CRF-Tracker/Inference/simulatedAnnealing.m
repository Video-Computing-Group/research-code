function [stateFound] = simulatedAnnealing(adj,J,local,temperatures,schedule)
% [stateFound] = simulatedAnnealing(adj,J,local,temperatures,schedule)
% samples between states in the system, while cooling the system in the
% given cooling schedule.
%
% Input Arguments:
% 1. adj: NxN matrix, defines neighbours in the graph
% 2. J: NxN matrix, defines the strength of the interactions between the nodes.
%       if Psi{i,j} is the general MRF matrix with the potentials between
%       i and j, than:
%       Psi{i,j} = [exp(J(i,j)), 1; 1, exp(J(i,j))];
% 3. local: 1xN cell array, each cell is a column-vector, Vix1, defines the
%           local potentials for each node.
% 4. temperatures, schedule: both row vectors of the same length, representing
%                            the cooling schedule. schedule(i) is the number of
%                            steps to be taken at temperatures(i)
%
% Output Arguments:
% stateFound: the state the system reached at the end of the sampling
%             process
%

adjCell = adjMat2Cell(adj);
N = length(local);
lambda = cell(1,N);
for i=1:N
    lambda{i} = J(i,adjCell{i});
end
stateFound = sim_anneal(adjCell,lambda,local,temperatures,schedule);
