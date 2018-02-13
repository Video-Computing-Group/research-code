% Implementation of a (directed) cut function
% Author: Andreas Krause (krausea@gmail.com)
%
% function C = sfo_fn_cutfun(G,A)
% G: Adjacency matrix of the graph
% A: subset of vertices to measure cut *from*
%
% Example: G = [1 1 0; 1 0 1; 0 1 1]; F = sfo_fn_cutfun(G); F([1 3])

function F = sfo_fn_cutfun(G)
fn = @(A) eval_cut_fn(G,A);
F = sfo_fn_wrapper(fn);


function C = eval_cut_fn(G,A)
A = sfo_unique_fast(A);
n = size(G,1);
C = sum(sum(G(A,sfo_setdiff_fast(1:n,A))));
