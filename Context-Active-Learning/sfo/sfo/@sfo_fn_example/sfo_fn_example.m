% The example from the tutorial slides at www.submodularity.org
% Implemented by Andreas Krause (krausea@gmail.com)
% F([]) = 0, F([1])= -1, F([2]) = 2, F([1,2]) = 0
%
% function R = sfo_fn_example(A)
% A: Input set to evaluate, [],[1],[2],[1,2]
% Example: F = sfo_fn_example; F([1,2])

function F = sfo_fn_example
F = sfo_fn_wrapper(@example_fun);

function R = example_fun(A)

load SFMinData
x = zeros(size(SFMinData.M,1),1);
x(A) = 1;
% R = sum(SFMinData.M*x) + SFMinData.gamma*length(A);
R = 0.5*x'*SFMinData.Q*x + x'*SFMinData.f + SFMinData.gamma*length(A);