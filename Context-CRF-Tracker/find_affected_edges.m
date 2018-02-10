function [affected_edges] = find_affected_edges(consolidated_el, ednum, correspondence_mat)

affected_edges = [];
t1 = consolidated_el(ednum,1);
t2 = consolidated_el(ednum,2);
z1 = consolidated_el(ednum,3);
z2 = consolidated_el(ednum,4);
n1 = consolidated_el(ednum,5);
n2 = consolidated_el(ednum,6);

cmat = correspondence_mat;
cmat{t1,t2}{z1,z2}(n1,n2) = 1;

if [t1 t2 z1 z2] == [1 1 1 2]
    cm = cmat{1,1}{1,2}*cmat{1,2}{2,2}*cmat{2,2}{1,2}';
    [I,J] = find(cm > correspondence_mat{1,2}{1,1});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 2 1 1 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,2}{1,1}'*cmat{1,1}{1,2}*cmat{1,2}{2,2};
    [I,J] = find(cm > correspondence_mat{2,2}{1,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [2 2 1 2 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,1}{1,2}'*cmat{1,2}{1,1}*cmat{2,2}{1,2};
    [I,J] = find(cm > correspondence_mat{1,2}{2,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 2 2 2 I(i) J(i)]];
        end;
        return;
    end;
elseif [t1 t2 z1 z2] == [1 2 1 1]
    cm = cmat{1,2}{1,1}*cmat{2,2}{1,2}*cmat{1,2}{2,2}';
    [I,J] = find(cm > correspondence_mat{1,1}{1,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 1 1 2 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,2}{1,1}'*cmat{1,1}{1,2}*cmat{1,2}{2,2};
    [I,J] = find(cm > correspondence_mat{2,2}{1,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [2 2 1 2 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,1}{1,2}'*cmat{1,2}{1,1}*cmat{2,2}{1,2};
    [I,J] = find(cm > correspondence_mat{1,2}{2,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 2 2 2 I(i) J(i)]];
        end;
        return;
    end;
elseif [t1 t2 z1 z2] == [2 2 1 2]
    cm = cmat{1,2}{1,1}*cmat{2,2}{1,2}*cmat{1,2}{2,2}';
    [I,J] = find(cm > correspondence_mat{1,1}{1,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 1 1 2 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,1}{1,2}*cmat{1,2}{2,2}*cmat{2,2}{1,2}';
    [I,J] = find(cm > correspondence_mat{1,2}{1,1});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 2 1 1 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,1}{1,2}'*cmat{1,2}{1,1}*cmat{2,2}{1,2};
    [I,J] = find(cm > correspondence_mat{1,2}{2,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 2 2 2 I(i) J(i)]];
        end;
        return;
    end;
elseif [t1 t2 z1 z2] == [1 2 2 2]
    cm = cmat{1,2}{1,1}*cmat{2,2}{1,2}*cmat{1,2}{2,2}';
    [I,J] = find(cm > correspondence_mat{1,1}{1,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 1 1 2 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,1}{1,2}*cmat{1,2}{2,2}*cmat{2,2}{1,2}';
    [I,J] = find(cm > correspondence_mat{1,2}{1,1});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [1 2 1 1 I(i) J(i)]];
        end;
        return;
    end;
    cm = cmat{1,2}{1,1}'*cmat{1,1}{1,2}*cmat{1,2}{2,2};
    [I,J] = find(cm > correspondence_mat{2,2}{1,2});
    if ~isempty(I) & ~isempty(J)
        for i = 1:length(I)
            affected_edges = [affected_edges; [2 2 1 2 I(i) J(i)]];
        end;
        return;
    end;
end;


    

