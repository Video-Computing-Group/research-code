function prbmat = prob_mat(Dmat,thresh)
%%%%% thresh can be vector or a scalar value
 Dm=Dmat;
if (length(thresh)==1)
  temp0=1./(1+exp(Dm-thresh)) ; 
  prbmat=temp0;  
else
   temp0=repmat(thresh,1,size(Dm,2));
   temp1=1./(1+exp(Dm-temp0)) ; 
   prbmat=temp1;
     
end







end