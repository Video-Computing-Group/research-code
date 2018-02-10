function D = mhd_computation_pt_cloud(cloud, A, x)
D = Inf.*ones(size(A,1), size(cloud,1));
progressbar('Creating MHD Map');
for c = 1:size(A,1)
    if ~isempty(x{c}) & ~isempty(A{c})
        D(c,:) = sum(((cloud - ones(length(cloud),1)*x{c}')*A{c}).*(cloud - ones(length(cloud),1)*x{c}'),2)';
        progressbar(c/size(A,1));
    end;
end;
return;