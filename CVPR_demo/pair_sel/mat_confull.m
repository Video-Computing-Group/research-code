function [A_ineqportion]= mat_confull(num_p)
num=num_p^2;
camNamearr=[1 2 3 4];
triplet=nchoosek(camNamearr,3);
pairs=nchoosek(camNamearr,2);
pairs_merged=str2num(strcat(num2str(pairs(:,1)),num2str(pairs(:,2))));


load con_mat.mat
con_mat=sparse(con_mat);
%init_mat=zeros(size(con_mat,1)*size(triplet,1),size(pairs,1));
init_mat2=[];
part_12=con_mat(:,1:num);
part_13=con_mat(:,num+1:2*num);
part_23=con_mat(:,2*num+1:end);
part_zeros=zeros(size(part_23));
for i=1:size(triplet,1)
init_mat1=[];
ct_pairs=nchoosek(triplet(i,:),2);
ct_pairs_merged=str2num(strcat(num2str(ct_pairs(:,1)),num2str(ct_pairs(:,2))));
[loct,~]=find(pairs_merged==ct_pairs_merged');
count=1;
    for j=1:size(pairs,1)
        
        
        
        
        if(((sum(loct==j))~=0) && count==1)
            temp_mat=part_12;
            count=count+1;
        elseif ((sum(loct==j))~=0 && count==2)
            temp_mat=part_13;
            count=count+1;

        elseif((sum(loct==j))~=0 && count==3)
            temp_mat=part_23;
            count=1;
        end
        
        if((sum(loct==j))==0 )
            temp_mat=part_zeros;
             %count=1;
        end
       init_mat1=[init_mat1 temp_mat]; 
        
    end
init_mat2=[init_mat2;init_mat1];    
    
end

A_ineqportion=sparse(init_mat2);




%%



end