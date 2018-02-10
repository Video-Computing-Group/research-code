function [ep] = compute_edge_potential(cord_1, cord_2, cell_1, cell_2, states_1, states_2)

k = length(states_1);
ep = zeros(k+1);
ep(1,:) = 1/((k+1)^2).*ones(1, (k+1));
ep(2:end, 1) = 1/((k+1)^2).*ones(k, 1);
rest = 1 - sum(ep(1,:)) - sum(ep(2:end, 1));

c_1 = mean(fillInside(cord_1(cell_1).PixelList), 1);
c_2 = mean(fillInside(cord_1(cell_2).PixelList), 1);
c_1_c_2 = c_2 - c_1;

row = 2;
for i = states_1
    col = 2;
    s_1 = mean(fillInside(cord_2(i).PixelList), 1);
    for j = states_2
        s_2 = mean(fillInside(cord_2(j).PixelList), 1);
        s_1_s_2 = s_2 - s_1;
        ep(row, col) = exp(-norm(s_1_s_2 - c_1_c_2));
        col = col + 1;
    end;
    row = row + 1;
end;

ep(2:end, 2:end) = ep(2:end, 2:end) .* (rest/sum(sum(ep(2:end, 2:end))));

return;
        
