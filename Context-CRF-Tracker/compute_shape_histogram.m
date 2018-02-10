function [N] = compute_shape_histogram(points, method, min_r, max_r, min_th, max_th, n_r, n_th)

[th, r] = cart2pol(points(:,1), points(:,2));

if strcmpi(method, 'shell')
    r_edges = [min_r : (max_r-min_r)/n_r : max_r];
    N  = histc(r, r_edges);
    N = N./sum(N);
elseif strcmpi(method, 'weighted_shell')
    N = zeros(1, n_r);
    count = 1;
    for i = min_r:(max_r-min_r)/n_r:max_r-(max_r-min_r)/n_r
        N(count) = mean(th(find(r>=i & r<i+(max_r-min_r)/n_r)));
        count = count + 1;
    end;
elseif strcmpi(method, 'normalized_weighted_shell')
    N = zeros(1, n_r);
    count = 1;
    for i = min_r:(max_r-min_r)/n_r:max_r-(max_r-min_r)/n_r
        N(count) = mean(th(find(r>=i & r<i+(max_r-min_r)/n_r)));
        count = count + 1;
    end;
    N = N./sum(N);
    
elseif strcmpi(method, 'sector')
    th_edges = [min_th : (max_th-min_th)/n_th : max_th];
    N  = histc(th, th_edges);
    N = N./sum(N);
elseif strcmpi(method, 'weighted_sector')
    N = zeros(1, n_th);
    count = 1;
    for i = min_th:(max_th-min_th)/n_th:max_th-(max_th-min_th)/n_th
        N(count) = mean(r(find(th>=i & th<i+(max_th-min_th)/n_th)));
        count = count + 1;
    end;
elseif strcmpi(method, 'normalized_weighted_sector')
    N = zeros(1, n_th);
    count = 1;
    for i = min_th:(max_th-min_th)/n_th:max_th-(max_th-min_th)/n_th
        N(count) = mean(r(find(th>=i & th<i+(max_th-min_th)/n_th)));
        count = count + 1;
    end;
    N = N./sum(N);
    
elseif strcmpi(method, 'spiderweb')
    N = zeros(1, n_th*n_r);
    count = 1;
    for i = min_r:(max_r-min_r)/n_r:max_r-(max_r-min_r)/n_r
        for j = min_th:(max_th-min_th)/n_th:max_th-(max_th-min_th)/n_th
            N(count) = sum(r>=i & r<i+(max_r-min_r)/n_r & th>=j & th<j+(max_th-min_th)/n_th);
            count = count + 1;
        end;
    end;
    N = N./sum(N);
end;

return;

