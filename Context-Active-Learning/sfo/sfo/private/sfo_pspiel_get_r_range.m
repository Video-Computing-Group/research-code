% Helper function for sfo_pspiel, by Andreas Krause (krausea@gmail.com)
% Example: See sfo_tutorial.m
%% get a reasonable range of Rs from shortest path matrix
function Rs = sfo_pspiel_get_r_range(spdist,k)
k = k-1;
rng = sort(spdist(:));
N = length(rng);
Rs = zeros(1,k);
for i=1:k
    ind = floor((N-floor(sqrt((k-i)/k)*N))*.75)+1;
    Rs(i)= rng(ind);
end
% ensure last value of Rs leads to "trivial" padded decomposition
Rs = [Rs 2*max(spdist(:))];

