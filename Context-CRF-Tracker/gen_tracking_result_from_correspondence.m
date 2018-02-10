function [tracking_result, untracked, division_result] = gen_tracking_result_from_correspondence(correspondence)

tracking_result = [];
untracked = [];
division_result = [];

for i = 1 : size(correspondence, 1)
    if isempty(find(correspondence(i,:)==1))
        untracked = [untracked; i];
    elseif length(find(correspondence(i,:)==1)) == 1
        tracking_result = [tracking_result; [i, find(correspondence(i,:)==1)]];
    else
        division_result = [division_result; [i, find(correspondence(i,:)==1)]];
    end;
end;

return;