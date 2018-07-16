function ineqportion1 = pwise_wbasedsel(np_array,nv)
pairs=nchoosek((1:length(np_array)),2);
ineqportion1=zeros(size(pairs,1),nv);
bias=0;
for i=1:length(pairs)
  ct_cama=pairs(i,1);
  ct_camb=pairs(i,2);
  ct_pcount= np_array(ct_cama)*np_array(ct_camb);
  ineqportion1(i,bias+1:ct_pcount+bias)=1;
  bias=bias+ct_pcount;
  
    
    
    
    
end







end