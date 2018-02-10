function [assoc] = GNCDA2GRP(pi, kk)

tic;
cameras = [1, 2];
numPersons = zeros(1, length(cameras));
for c = 1 : size(pi,1) - 1
    numPersons(c) = size(pi{c,c+1},1);
end;
numPersons(size(pi,1)) = size(pi{size(pi,1)-1,size(pi,1)},2);

CPairs = nchoosek(cameras,2);
numCPairs = size(nchoosek(cameras,2), 1);
numCameras = length(cameras);

model = struct('obj', [], 'A', sparse([]), 'rhs', [], 'sense', [], 'modelsense', 'min', 'vtype', 'B');

vectorSize_perpair = zeros(numCPairs, 1);
for i = 1:numCPairs
    n1 = CPairs(i,1);
    n2 = CPairs(i,2);
    vectorSize_perpair(i) = numPersons(n1)*numPersons(n2);
    cc =  pi{n1,n2}(1:numPersons(n1),1:numPersons(n2));
    model.obj = [model.obj; cc(:)];
end;
model.obj = kk.*ones(length(model.obj),1) - model.obj;
disp('Cost vector done.');

triplets = [];
vectorSize_pertriplet = [];
vectorSize_perquartlet = [];

vectorSize = size(model.obj,1);
m = (numCameras-1)*sum(numPersons) + sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet);
n = vectorSize;
I = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)));
J = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)));
S = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)));
model.rhs = ones(m, 1);
counter = 1;
rowcount = 1;
for c = 1:numCPairs
    n1 = CPairs(c,1);
    n2 = CPairs(c,2);    
    for i = 1:numPersons(n1)
        for j = 1:numPersons(n2)
            I(1, counter) = rowcount;
            J(1, counter) = sum(vectorSize_perpair(1:c-1))+(j-1)*numPersons(n1)+i;
            S(1, counter) = 1;
            counter = counter + 1;
        end;
        rowcount = rowcount + 1;
    end;
    for j = 1:numPersons(n2)
        for i = 1:numPersons(n1)
            I(1, counter) = rowcount;
            J(1, counter) = sum(vectorSize_perpair(1:c-1))+(j-1)*numPersons(n1)+i;
            S(1, counter) = 1;
            counter = counter + 1;
        end;
        rowcount = rowcount + 1;
    end;
end;
disp('Pairwise constraint matrix done.');

model.A = sparse(I, J, S, m, n);
clear I J S m n;
model.sense = '<';
model.modelsense = 'min';
model.vtype = 'B';

%%
result = cplexbilp(model.obj, model.A, model.rhs, [], []);
clear model

%%
assoc = cell(numCameras,numCameras);
for cp = 1:numCPairs
    ci = CPairs(cp,1);
    cj = CPairs(cp,2);
    assoc{ci,cj} = reshape(result(sum(vectorSize_perpair(1:cp-1))+1:sum(vectorSize_perpair(1:cp))),[numPersons(ci) numPersons(cj)]);
end;

return;