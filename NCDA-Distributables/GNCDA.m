function [assoc] = GNCDA(pi, kk)

cameras = [1 : size(pi,1)];
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
for r = 1:size(CPairs,1)
    ci = CPairs(r,1);
    ck = CPairs(r,2);
    others = setdiff(cameras,[ci,ck]);
    all_perms_tri = nchoosek(others, 1);
    for cj = all_perms_tri'
        triplets = [triplets; [ci,cj,ck]];
        vectorSize_pertriplet = [vectorSize_pertriplet; ...
            [numPersons(ci)*numPersons(ck)*numPersons(cj)]];
    end;
end;

vectorSize = size(model.obj,1);
m = (numCameras-1)*sum(numPersons) + sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet);
n = vectorSize;
I = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)));
J = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)));
S = zeros(1, sum(2.*vectorSize_perpair) + 3*(sum(vectorSize_pertriplet) + sum(vectorSize_perquartlet)));
model.rhs = 2.*ones(m, 1);
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

for t = 1:size(triplets,1)
    ci = triplets(t,1);
    cj = triplets(t,2);
    ck = triplets(t,3);
    cij = find((CPairs(:,1)==ci & CPairs(:,2)==cj) | (CPairs(:,1)==cj & CPairs(:,2)==ci));
    cjk = find((CPairs(:,1)==cj & CPairs(:,2)==ck) | (CPairs(:,1)==ck & CPairs(:,2)==cj));
    cik = find((CPairs(:,1)==ci & CPairs(:,2)==ck) | (CPairs(:,1)==ck & CPairs(:,2)==ci));
    for i = 1:numPersons(ci)
        for k = 1:numPersons(ck)
            for j = 1:numPersons(cj)
                if ci < cj
                    I(1, counter) = rowcount;
                    J(1, counter) = sum(vectorSize_perpair(1:cij-1))+(j-1)*numPersons(ci)+i;
                    S(1, counter) = 1;
                    counter = counter + 1;
                else
                    I(1, counter) = rowcount;
                    J(1, counter) = sum(vectorSize_perpair(1:cij-1))+(i-1)*numPersons(cj)+j;
                    S(1, counter) = 1;
                    counter = counter + 1;
                end;
                if cj < ck
                    I(1, counter) = rowcount;
                    J(1, counter) = sum(vectorSize_perpair(1:cjk-1))+(k-1)*numPersons(cj)+j;
                    S(1, counter) = 1;
                    counter = counter + 1;
                else
                    I(1, counter) = rowcount;
                    J(1, counter) = sum(vectorSize_perpair(1:cjk-1))+(j-1)*numPersons(ck)+k;
                    S(1, counter) = 1;
                    counter = counter + 1;
                end;
                if ci < ck
                    I(1, counter) = rowcount;
                    J(1, counter) = sum(vectorSize_perpair(1:cik-1))+(k-1)*numPersons(ci)+i;
                    S(1, counter) = -1;
                    counter = counter + 1;
                else
                    I(1, counter) = rowcount;
                    J(1, counter) = sum(vectorSize_perpair(1:cik-1))+(i-1)*numPersons(ck)+k;
                    S(1, counter) = -1;
                    counter = counter + 1;
                end;
                rowcount = rowcount + 1;
            end;
        end;
    end;
end;
model.A = sparse(I, J, S, m, n);
clear I J S m n;
disp('Loop constraint matrix done.');
model.sense = '<';
model.modelsense = 'min';
model.vtype = 'B';

%% Solve NCDA
result = cplexbilp(model.obj, model.A, model.rhs, [], []);
clear model

%% Generate final association matrix from output of optimization (vector)
assoc = cell(numCameras,numCameras);
for cp = 1:numCPairs
    ci = CPairs(cp,1);
    cj = CPairs(cp,2);
    assoc{ci,cj} = reshape(result(sum(vectorSize_perpair(1:cp-1))+1:sum(vectorSize_perpair(1:cp))),[numPersons(ci) numPersons(cj)]);
end;

return;