function [candidates] = generate_division_candidates(cord_1, cord_2, adj, properties)
%#codegen
%coder.inline('Never')

candidates = [];


for c1 = 1 : length(cord_1)
    nbors = get_k_closest_neighbors(cord_1, c1, cord_2, properties.Number_of_Nonzero_States + 1);
    pairs = combntns(nbors, 2);
    for i = 1 : double(size(pairs, 1))
        if adj(pairs(i,1), pairs(i,2)) == 1
            [isdivided, d1, d2] = detect_cell_divisions(cord_1, c1, cord_2, pairs(i,:), ...
                [properties.Division_MHD_Threshold properties.Division_Area_Threshold], 1, 1);
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

return;