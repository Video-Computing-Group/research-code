function rho = rhoMat(rhoVal,adj)

rhoSize = size(rhoVal);
N = size(adj,1);
neighbNum = sum(adj);

if (iscell(rhoVal))
    error('argument given for trw (in loopy) must either be a scalar or a matrix in size <num_nodes> x <num_nodes>');
end

if prod(rhoSize)==1
    rho = cell(1,N);
    for i=1:N
        rho{i} = rhoVal*ones(1,neighbNum(i));
    end
else
    if (length(rhoSize)~=2 || rhoSize(1)~=N || rhoSize(2)~=N)
        error('argument given for trw (in loopy) must either be a scalar or a matrix in size <num_nodes> x <num_nodes>');
    else
        rho = cell(1,N);
        for i=1:N
            if neighbNum(i)>0
                rho{i} = rhoVal(i,find(adj(i,:)));
            end
        end
    end
end
