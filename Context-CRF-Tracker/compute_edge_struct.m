function [edge_struct] = compute_edge_struct(adj, states)

nStates = [];
for i = 1:length(states)
    nStates(i,1) = int32(length(states(i,:))+1);
end;

edge_struct = UGM_makeEdgeStruct(int32(adj),nStates);

return;