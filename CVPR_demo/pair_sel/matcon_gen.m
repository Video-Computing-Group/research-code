function out= matcon_gen(m,n,p)

count=1;
for i=1:m*n
 temp_ab=zeros(m*n,1);
 temp_ab(i,1)=1;
 temp_m=reshape(temp_ab,m,n);
 
 for j=1:p
     [row_a,col_a]=find(temp_m==1);
     temp_ac=zeros(m,p);
     temp_bc=zeros(n,p);
     temp_ac(row_a,j)=1;
     temp_bc(col_a,j)=1;
     out(count,:)= [(temp_ab)' (temp_ac(:))' (temp_bc(:))'];
     count=count+1;
     
 end
 
 
    
    
    
    
    
    
    
end







end