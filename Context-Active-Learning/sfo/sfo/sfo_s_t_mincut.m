% Finding the minimum A of a submodular function such that s in A and t not in A
% Implementation by Andreas Krause (krausea@gmail.com)
%
% function A = sfo_s_t_mincut(F,V,s,t)
% F: Submodular function
% V: index set
% s: element in V to include
% t: element in V to exclude
% opt (optional): options, depending on 
% 
% minnorm_stopping_thresh: threshold for stopping minimization (1e-5)
%
% Example: 
%   G = [1 1 0; 1 0 1; 0 1 1]; F = sfo_fn_cutfun(G)
%   A = sfo_s_t_mincut(F,1:3,2,1)


function A = sfo_s_t_mincut(F,V,s,t,opt)
if ~exist('opt','var')
    opt = sfo_opt({'minnorm_stopping_thresh',1e-5});
end
% F2 guarantees we pick s
F2 = sfo_fn_residual(F,s);
% V2 guarantees we can't pick s and t
V2 = sfo_setdiff_fast(V,[s,t]);
A = sfo_min_norm_point(F2,V2,opt);
A = [A s];
