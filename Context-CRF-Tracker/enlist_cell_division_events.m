function [division_results] = enlist_cell_division_events(cord_1, cord_2, properties) 
%#codegen
division_results = [];

[adj] = compute_adjacency_matrix(cord_2);
candidates = [];

for c1 = 1 : length(cord_1)
    nbors = get_k_closest_neighbors(cord_1, c1, cord_2, properties.Number_of_Nonzero_States + 1);
    pairs = combntns(nbors, 2);
    for i = 1 : size(pairs, 1)
        if adj(pairs(i,1), pairs(i,2)) == 1
            [isdivided, d1, d2] = detect_cell_divisions(cord_1, c1, cord_2, pairs(i,:), ...
                [properties.Division_MHD_Threshold properties.Division_Area_Threshold]);
            if isdivided
                candidates = [candidates; [c1 pairs(i,:) ...
                    (d1/properties.Division_MHD_Threshold + d2/properties.Division_Area_Threshold)]];
            end;
        end;
    end;
end;
       
if size(candidates, 1) >= 2
    candidates = sortrows(candidates, 4);
end;

while size(candidates, 1) > 0
    division_results = [division_results; candidates(1, 1:3)];
    I1 = find(candidates(:,1) == candidates(1,1));
    I2 = find(candidates(:,2) == candidates(1,2));
    I3 = find(candidates(:,3) == candidates(1,2));
    I4 = find(candidates(:,2) == candidates(1,3));
    I5 = find(candidates(:,3) == candidates(1,3));
    I = unique([I1; I2; I3; I4; I5]);
    candidates(I,:) = [];
end;

return;