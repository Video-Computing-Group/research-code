function vertex_set =vertex_set_gen(np_array)
vertex_set=cell(1,length(np_array));

vertex_set{1}=1:np_array(1);
 for i=2:length(np_array)
    temp1=1:np_array(i);
    temp2= vertex_set{i-1}(end);
    vertex_set{i}=temp1+temp2;
    
    
    
 end





end