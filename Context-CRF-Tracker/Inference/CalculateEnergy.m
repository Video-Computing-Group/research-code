function E = CalculateEnergy(gridM, lambdaM, localM, assign)

epsilon = 1e-200;

% E = 1;
% nNodes = length(gridM);
% for j=1:nNodes, 
%     E = E * localM{j}(assign(j));
%     for j1=j+1:nNodes,
%         if(gridM(j,j1))
%             E = E * lambdaM{j,j1}(assign(j),assign(j1));
%         end
%     end
% end
% 
% E = -E;
% 
E = 0;
nNodes = length(gridM);
for j=1:nNodes,
    E = E - log(localM{j}(assign(j)) + epsilon);
    N = find(gridM(j,:));
    for j1=1:length(N),
        if (N(j1)<j) continue; end
        E = E - log(lambdaM{j,N(j1)}(assign(j), assign(N(j1)))+epsilon);
    end
end