function [correspondence_mat] = cell_tracker(cord, mode, properties_spat, properties_temp)

if strcmpi(mode, 'spatial')
    
    correspondence_mat = cell(1,length(cord));
    for i = 1 : length(cord)-1
        [c, el] = pairwise_similarity_score_generator(cord{i}, cord{i+1}, mode, properties_spat);
        total_cost = 0;
        el(:,3) = ones(size(el,1),1) - 2.*el(:,3);
        el = sortrows(el, [3 -1]);
        while ~isempty(el) & el(1,3) < 0
            n1 = el(1,1);
            n2 = el(1,2);
            total_cost = total_cost + el(1,3);
            c(n1,n2) = 1;
            el = remove_edges_from_list(el, n1, n2);
        end;
        correspondence_mat{i+1} = c;
    end;

elseif strcmpi(mode, 'temporal')
    
    correspondence_mat = cell(1,length(cord));
    for i = 1 : length(cord)-1
        [c, el] = pairwise_similarity_score_generator(cord{i}, cord{i+1}, mode, properties_temp);
        total_cost = 0;
        el(:,3) = ones(size(el,1),1) - 2.*el(:,3);
        el = sortrows(el, [3 -1]);
        while ~isempty(el) & el(1,3) < 0
            n1 = el(1,1);
            n2 = el(1,2);
            total_cost = total_cost + el(1,3);
            c(n1,n2) = 1;
            el = remove_edges_from_list(el, n1, n2);
        end;
        correspondence_mat{i+1} = c;
    end;
    
elseif strcmpi(mode, 'spatio-temporal')
    
    correspondence_mat = cell(length(cord), length(cord));
    correspondence_mat{1,1} = cell(length(cord{1}), length(cord{1}));
    correspondence_mat{2,2} = cell(length(cord{2}), length(cord{2}));
    correspondence_mat{1,2} = cell(length(cord{1}), length(cord{2}));
    el = cell(length(cord), length(cord));
    el{1,1} = cell(length(cord{1}), length(cord{1}));
    el{2,2} = cell(length(cord{2}), length(cord{2}));
    el{1,2} = cell(length(cord{1}), length(cord{2}));
    list_of_cell_division_edges = [];
    for i = 1 : length(cord{1})-1
        [correspondence_mat{1,1}{i,i+1}, el{1,1}{i,i+1}] = pairwise_similarity_score_generator(cord{1}{i}, cord{1}{i+1}, 'spatial', properties_spat);
        el{1,1}{i,i+1}(:,3) = ones(size(el{1,1}{i,i+1},1),1) - 2.*el{1,1}{i,i+1}(:,3);
    end;
    for i = 1 : length(cord{2})-1
        [correspondence_mat{2,2}{i,i+1}, el{2,2}{i,i+1}] = pairwise_similarity_score_generator(cord{2}{i}, cord{2}{i+1}, 'spatial', properties_spat);
        el{2,2}{i,i+1}(:,3) = ones(size(el{2,2}{i,i+1},1),1) - 2.*el{2,2}{i,i+1}(:,3);
    end; 
    for i = 1 : min(length(cord{1}),length(cord{2}))
        [correspondence_mat{1,2}{i,i}, el{1,2}{i,i}, div_list] = pairwise_similarity_score_generator(cord{1}{i}, cord{2}{i}, 'temporal', properties_temp); % Modified
        list_of_cell_division_edges = [list_of_cell_division_edges; [ones(size(div_list,1),1)*[1 2 i i], div_list]]; % New
        el{1,2}{i,i}(:,3) = ones(size(el{1,2}{i,i},1),1) - 2.*el{1,2}{i,i}(:,3);
    end;
    
    % Similarity socres generated and division edges set aside

    
    c = cell(2,2);
    c{1,1} = cell(2,2);
    c{2,2} = cell(2,2);
    c{1,2} = cell(2,2);
    e = cell(2,2);
    e{1,1} = cell(2,2);
    e{2,2} = cell(2,2);
    e{1,2} = cell(2,2);
    c{1,1}{1,2} = correspondence_mat{1,1}{1,2};
    c{2,2}{1,2} = correspondence_mat{2,2}{1,2};
    c{1,2}{1,1} = correspondence_mat{1,2}{1,1};
    c{1,2}{2,2} = correspondence_mat{1,2}{2,2};
    e{1,1}{1,2} = el{1,1}{1,2};
    e{2,2}{1,2} = el{2,2}{1,2};
    e{1,2}{1,1} = el{1,2}{1,1};
    e{1,2}{2,2} = el{1,2}{2,2};
    c = data_association_atomic_quartet(e, c);
    correspondence_mat{1,1}{1,2} = c{1,1}{1,2};
    correspondence_mat{2,2}{1,2} = c{2,2}{1,2};
    correspondence_mat{1,2}{1,1} = c{1,2}{1,1};
    correspondence_mat{1,2}{2,2} = c{1,2}{2,2};
    
    for z = 2 : length(cord{1})-1
        c = cell(2,2);
        c{1,1} = cell(2,2);
        c{2,2} = cell(2,2);
        c{1,2} = cell(2,2);
        e = cell(2,2);
        e{1,1} = cell(2,2);
        e{2,2} = cell(2,2);
        e{1,2} = cell(2,2);
        c{1,1}{1,2} = correspondence_mat{1,1}{z,z+1};
        c{2,2}{1,2} = correspondence_mat{2,2}{z,z+1};
        c{1,2}{1,1} = correspondence_mat{1,2}{z,z};
        c{1,2}{2,2} = correspondence_mat{1,2}{z+1,z+1};
        e{1,1}{1,2} = el{1,1}{z,z+1};
        e{2,2}{1,2} = el{2,2}{z,z+1};
        e{1,2}{2,2} = el{1,2}{z+1,z+1};
        c = data_association_atomic_quartet(e, c);
        correspondence_mat{1,1}{z,z+1} = c{1,1}{1,2};
        correspondence_mat{2,2}{z,z+1} = c{2,2}{1,2};
        correspondence_mat{1,2}{z,z} = c{1,2}{1,1};
        correspondence_mat{1,2}{z+1,z+1} = c{1,2}{2,2};
    end;
    
     % Division edges re-incorporated

    
    for i = 1 : size(list_of_cell_division_edges,1)
        t1 = list_of_cell_division_edges(i,1);
        t2 = list_of_cell_division_edges(i,2);
        z1 = list_of_cell_division_edges(i,3);
        z2 = list_of_cell_division_edges(i,4);
        n1 = list_of_cell_division_edges(i,5);
        n2 = list_of_cell_division_edges(i,6);
        correspondence_mat{t1,t2}{z1,z2}(n1,n2)=1;
    end;
    
end;
    
return;


