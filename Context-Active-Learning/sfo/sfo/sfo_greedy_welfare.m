% Andreas Krause (krausea@gmail.com)
% greedy algorithm for submodular welfare 
% (maximization on partition matroid)
% Solves the problem
% [A1',...,Ak'] = max_{A1,...,Ak} sum_i F_i(A_i)
%
% function A_part = sfo_greedy_welfare(Fs,V,k)
% Fs: cell array of monotonic submodular functions
% V: index set
% k: #elements to pick (optional)
%
% Example: (also see sfo_tutorial.m)
%   F = sfo_fn_varred(0.5*eye(10)+0.5*ones(10),1:10)
%   Fs = {F,F}
%   A_part = sfo_greedy_welfare(Fs,1:10)


function [A_part,scores] = sfo_greedy_welfare(Fs,V,k)
if ~exist('k','var')
    k = length(V);
end
m = length(Fs); %number of buckets
n = length(V); % number of elements
% we generate a new set that contains "colored" versions of the original
% elements. with this construction, for each element i in V, mod(i,m) is
% the bucket id (color), and floor(i/m) is the original item number
V_welfare = repmat([0:(m-1)]',1,n)+ones(m,1)*V*m; V_welfare = V_welfare(:)';
% construct the utility function that deals with the colored elements

F_welfare = sfo_fn_welfare(Fs);


% need to specify a function that identifies the independent sets of the
% partition matroid
opt = sfo_opt({'greedy_check_indep',(@(As) sfo_check_welfare(Fs,As))});

% now we just run the basic greedy algorithm
A_welfare = sfo_greedy_lazy(F_welfare,V_welfare,k,opt);

% and "reinterpret" the solution by putting elements with the same "color"
% in a bucket
A_part = partition(F_welfare,A_welfare);

% return scores
scores = zeros(1,m);
for i = 1:m
    scores(i) = Fs{i}(A_part{i});
end

%% checks that no element is assigned two colors at the same time
function indep = sfo_check_welfare(Fs,As)
m = length(Fs);
els = floor(As/m);
indep = (length(els)==length(unique(els))); % a set is independent iff each el is picked at most once

