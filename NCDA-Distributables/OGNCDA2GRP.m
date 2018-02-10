function [assoc] = OGNCDA2GRP(pi, ass, kk)

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

CPairs_nonDummy = nchoosek(cameras(1:end-1), 2);
vectorSize_perpair_nonDummy = zeros(size(CPairs_nonDummy, 1), 1);
for i = 1:size(CPairs_nonDummy, 1)
    n1 = CPairs_nonDummy(i,1);
    n2 = CPairs_nonDummy(i,2);
    vectorSize_perpair_nonDummy(i) = numPersons(n1)*numPersons(n2);
end;

vectorSize = size(model.obj,1);
m = (numCameras-1)*sum(numPersons) + sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet) + sum(vectorSize_perpair_nonDummy); % The last sum(vectorSize_perpair_nonDummy)for the equality constraint  (already assigned labels)
n = vectorSize;
I = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)) + sum(vectorSize_perpair_nonDummy)); % The last sum(vectorSize_perpair_nonDummy)for the equality constraint 
J = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)) + sum(vectorSize_perpair_nonDummy)); % The last sum(vectorSize_perpair)for the equality constraint 
S = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)) + sum(vectorSize_perpair_nonDummy)); % The last sum(vectorSize_perpair)for the equality constraint 
model.rhs = ones(m, 1);
model.sense = blanks(m)';
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
        model.sense(rowcount, 1) = '<';
        rowcount = rowcount + 1;
    end;
    for j = 1:numPersons(n2)
        for i = 1:numPersons(n1)
            I(1, counter) = rowcount;
            J(1, counter) = sum(vectorSize_perpair(1:c-1))+(j-1)*numPersons(n1)+i;
            S(1, counter) = 1;
            counter = counter + 1;
        end;
        model.sense(rowcount, 1) = '<';
        rowcount = rowcount + 1;
    end;
end;
disp('Pairwise constraints done.');

for i = 1:size(CPairs_nonDummy, 1)
    ci = CPairs_nonDummy(i,1);
    cj = CPairs_nonDummy(i,2);
    cij = find((CPairs(:,1)==ci & CPairs(:,2)==cj) | (CPairs(:,1)==cj & CPairs(:,2)==ci));
    for i = 1 : numPersons(ci)
        for j = 1 : numPersons(cj)
            I(1, counter) = rowcount;
            J(1, counter) = sum(vectorSize_perpair(1:cij-1))+(j-1)*numPersons(ci)+i;
            S(1, counter) = 1;
            model.sense(rowcount, 1) = '=';
            model.rhs(rowcount, 1) = ass{ci,cj}(i,j);
            counter = counter + 1;
            rowcount = rowcount + 1;
        end;
    end;
end;
disp('Equality constraints done.');

model.A = sparse(I, J, S, m, n);
clear I J S m n;
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