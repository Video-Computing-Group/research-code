function [Obsclusters] =  gen_conncomp_4m_timestamp_xls(Timestamps_excel_file)                     %'Timestamps_6_Camera.xlsx'

timestamps = xlsread(Timestamps_excel_file);
numCameras = (size(timestamps,2)-1)/2;
counter = 1;
obslist = zeros(numCameras*size(timestamps, 1), 3);

for i = 1 : size(timestamps, 1)
    for cam = 1 : numCameras
        obslist(counter, 1) = 10*timestamps(i,1) + cam;
        obslist(counter, 2:3) = uint64(1000.*timestamps(i, 2*cam:2*cam+1));
        counter = counter + 1;
    end;
end;

%%
adj = adjacency_4m_timestamps(obslist);
[~,C] = graphconncomp(sparse(adj));
MC = cell(max(C), 1);
for c = 1 : max(C)
    MC{c} = find(C==c);
end;
Obsclusters = repmat(struct('Targets', [], 'TStart', [], 'TEnd', []), length(MC), 1);

for c = 1 : length(MC)
    obs = obslist(MC{c}',:);
    Obsclusters(c).Targets = repmat(struct('ID', [], 'Cam', [], 'TStart', [], 'TEnd', []), size(obs,1), 1);
    Obsclusters(c).TStart = min(obs(:,2));
    Obsclusters(c).TEnd = max(obs(:,3));
    for t = 1 : size(obs,1)
        Obsclusters(c).Targets(t).ID = floor(obs(t,1)/10);
        Obsclusters(c).Targets(t).Cam = mod(obs(t,1),10); 
        Obsclusters(c).Targets(t).TStart = obs(t,2); 
        Obsclusters(c).Targets(t).TEnd = obs(t,3);
    end;
end;

Obsclusters = sortStruct(Obsclusters, 'TStart');
Obsclusters = Obsclusters(2:end);

return;