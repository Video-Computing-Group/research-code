% Base class for set function objects
% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Functions are defined as objects representing set functions.
% For example, F = sfo_fn_entropy(sigma,1:size(sigma,1))
% will create a subclass of sfo_fn, such that F([1 4 6])
% will evaluate to the entropy of the Gaussian random variables
% with joint covariance sigma([1 4 6],[1 4 6]).
%
% The function objects are designed for efficiency. For example,
% computing the entropy of a set of D variables requires computing
% the Cholesky decomposition of a DxD matrix (O(n^3)). However,
% many algorithms for optimizing submodular functions require
% incremental computation, adding and removing elements to sets.
% The SFO function objects support this:
% Suppose A is a vector representing a set (such as [1 4 6]).
% Then FA = init(F,A) will return a new function object initialized
% with set A. This initialization has the advantage that it
% is very efficient to evaluate FA([A s]) for some new element s.
% This can be done using the inc method:
% val = inc(FA, A, s)
% will return the value FA([A s]).
% Note that even if FA is initialized with set A, the method
% inc(FA,B,s) will still correctly compute F([B s]), however
% it can not make use of efficiency gains by incremental evaluation.
%
% Set functions also support other operations, such as truncation,
% which are used in certain algorithms.
%
% Example: see the tutorial script for more detailed examples of the
% construction and use of set function objects.

function F = sfo_fn
F.current_set = -1;
F.current_val = 0;
F = class(F,'sfo_fn');
