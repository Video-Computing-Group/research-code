function [ep] = compute_edge_potential_gamma_temporal(cord_1, cord_2, cell_1, cell_2, states_1, states_2, gamma)

k = length(states_1);
ep = zeros(k);

c_1 = mean(fillInside(cord_1(cell_1).PixelList), 1);
c_2 = mean(fillInside(cord_1(cell_2).PixelList), 1);
c_1_c_2 = c_2 - c_1;
        
row = 1;
for i = states_1
    col = 1;
    s_1 = mean(fillInside(cord_2(i).PixelList), 1);
    for j = states_2
        s_2 = mean(fillInside(cord_2(j).PixelList), 1);
        s_1_s_2 = s_2 - s_1;
        ep(row, col) = exp(-gamma*norm(s_1_s_2 - c_1_c_2));
        col = col + 1;
    end;
    row = row + 1;
end;

return;