function [assoc, accuracy] = mainOGNCDA(Simscores, Obsclusters, kk)

numCameras = size(Simscores,1);
CPairs = nchoosek([1:numCameras], 2);

G = cell(numCameras,1);
assoc = cell(numCameras, numCameras);
GT = cell(numCameras, numCameras);
for i = 1 : size(assoc,1)
    for j = 1 : size(assoc,2)
        assoc{i,j} = zeros(size(Simscores{i,j},1), size(Simscores{i,j},2));
        GT{i,j} = zeros(size(Simscores{i,j},1), size(Simscores{i,j},2));
    end;
end;

accuracy = zeros(length(Obsclusters),1);

for o = 1 : length(Obsclusters)
    
    if o == 1
        for t = 1 : length(Obsclusters(o).Targets)
            tid = Obsclusters(o).Targets(t).ID;
            cam = Obsclusters(o).Targets(t).Cam;
            G{cam} = [G{cam}, tid];
        end;
        continue;
    end;
    
    D = [];
    CamD = [];
    for t = 1 : length(Obsclusters(o).Targets)
        tid = Obsclusters(o).Targets(t).ID;
        cam = Obsclusters(o).Targets(t).Cam;
        D = [D, tid];
        CamD = [CamD, cam];
    end;
    
    Fgrps = [];
    for cam = 1 : numCameras
        if ~isempty(G{cam})
            Fgrps = [Fgrps, cam];
        end;
    end;
    
    pi = cell(length(Fgrps)+1, length(Fgrps)+1);
    as = cell(length(Fgrps), length(Fgrps));
    
    for i = 1 : size(pi,1)-1
        for j = 1 : size(pi,1)-1
            if i >= j
                continue;
            end;
            cam1 = Fgrps(i);
            cam2 = Fgrps(j);
            pi{i,j} = zeros(length(G{cam1}), length(G{cam2}));
            as{i,j} = zeros(length(G{cam1}), length(G{cam2}));
            for m = 1 : length(G{cam1})
                for n = 1 : length(G{cam2})
                    tid1 = G{cam1}(m);
                    tid2 = G{cam2}(n);
                    pi{i,j}(m,n) = Simscores{cam1,cam2}(tid1,tid2);
                    as{i,j}(m,n) = assoc{cam1,cam2}(tid1,tid2);
                end;
            end;
        end;
    end;
    
    for i = 1 : size(pi,1)-1
        j = size(pi,1);
        cam1 = Fgrps(i);
        pi{i,j} = zeros(length(G{cam1}), length(D));
        for m = 1 : length(G{cam1})
            for n = 1 : length(D)
                cam2 = CamD(n);
                tid1 = G{cam1}(m);
                tid2 = D(n);
                if cam1 < cam2
                    pi{i,j}(m,n) = Simscores{cam1,cam2}(tid1,tid2);
                elseif cam1 > cam2
                    pi{i,j}(m,n) = Simscores{cam2,cam1}(tid2,tid1);
                else
                    pi{i,j}(m,n) = (-10^6);
                end;
            end;
        end;
    end;
    
    if length(Fgrps) < 2
        ass = GNCDA2GRP(pi, kk);
    else
        ass = OGNCDA(pi, as, kk);
    end;
    
    for i = 1 : size(ass,1)-1
        j = size(ass,1);
        cam1 = Fgrps(i);
        for m = 1 : length(G{cam1})
            for n = 1 : length(D)
                cam2 = CamD(n);
                tid1 = G{cam1}(m);
                tid2 = D(n);
                if cam1 < cam2
                    assoc{cam1,cam2}(tid1,tid2) = ass{i,j}(m,n);
                    GT{cam1,cam2}(tid1,tid2) = (tid1==tid2); %% Optional
                elseif cam1 > cam2
                    assoc{cam2,cam1}(tid2,tid1) = ass{i,j}(m,n);
                    GT{cam1,cam2}(tid1,tid2) = (tid1==tid2); %% Optional
                else
                    continue;
                end;
            end;
        end;
    end;
    
    for i = 1 : length(D)
        tid = D(i);
        cam = CamD(i);
        G{cam} = [G{cam}, tid];
    end;
    
    accc = 0;
    countt = 0;
    for l = 1 : size(CPairs,1)
        correct = 0;
        cam1 = CPairs(l,1);
        cam2 = CPairs(l,2);
        common = intersect(unique(G{cam1}), unique(G{cam2}));
        uncommon1 = setdiff(unique(G{cam1}), unique(G{cam2}));
        uncommon2 = setdiff(unique(G{cam2}), unique(G{cam1}));
        for i = common
            if assoc{cam1,cam2}(i,i) == 1
                correct = correct + 1;
            end;
        end;
        for i = uncommon1
            if sum(assoc{cam1,cam2}(i,:)) == 0
                correct = correct + 1;
            end;
        end;
        for i = uncommon2
            if sum(assoc{cam1,cam2}(:,i)) == 0
                correct = correct + 1;
            end;
        end;
        T = length(common) + length(uncommon1) + length(uncommon2);
        if T ~= 0
            accc = accc + correct/T;
            countt = countt + 1;
        end;
    end;
    accuracy(o) = accc/countt;
end;

return;
