function D = mhd_computation_pt_cloud_single_cell(cloud, A, x)
%#codegen
coder.inline('None');

D = Inf.*ones(1, size(cloud,1));
if ~isempty(x) && ~isempty(A)
    D = diag((cloud - ones(length(cloud),1)*x')*A*(cloud - ones(length(cloud),1)*x')')';
end;

return;