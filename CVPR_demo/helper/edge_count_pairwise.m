function out= edge_count_pairwise(pairs,np_array)
%%% gives pair/edge count per camera pair
for i=1:size(pairs,1)
    
    out(i)=np_array(pairs(i,1))*np_array(pairs(i,2));
    
end





end