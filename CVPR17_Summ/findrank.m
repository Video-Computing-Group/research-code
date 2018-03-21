% -------------------------------------------------------------------------
% Program to find the ranked list of shots by combining Z and Z_n
% Input: set of sparse coefficient matrices, Z and Z_n
% Output: ranked list of shots
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Ind = findrank(Z,t) 

N = size(Z,1);

r = zeros(1,N);
for i = 1:N
    r(i) = norm(Z(i,:),2);
end
[nrm,nrmInd] = sort(r,'descend');
nrmSum = 0;
for j = 1:N
    nrmSum = nrmSum + nrm(j);
    if (nrmSum / sum(nrm) > t)
        break;
    end
end
Ind = nrmInd(1:j);