% talyam
function i = chooseInteger(q)
% i = chooseInteger(q)
% given a row vector q (1XN), return i=1..N in probabilty q(i).
% for this purpose, I use rand to choose a number x, uniformly
% distributed on the interval (0.0, 1.0), and pick i in the
% following way:
% i = 1 if x<q(1)                      ==> Pr(1) = q(1)
% i = 2 if x>q(1) and x<(q(1)+q(2))    ==> Pr(2) = q(2)
%               ...                         ...
% i = N if x>(q(1)+...+q(N-1))         ==> Pr(N) = q(N)
x = rand;
N = length(q(1,:));
bins = cumsum(q);
i = 0;
for j=1:N
    if (x<bins(1,j))
        i = j;
        break;
    end
end
