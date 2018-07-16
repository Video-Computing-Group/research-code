function [A_ineqportion]= mat_confullg_bl(np_array)
camNamearr=1:length(np_array);
triplet=nchoosek(camNamearr,3);
pairs=nchoosek(camNamearr,2); %% [1-2 1-3 1-4 2-3 2-4 3-4]
pairs_merged=str2num(strcat(num_str_balance(pairs(:,1),length(np_array)),num_str_balance(pairs(:,2),length(np_array))));


%load con_mat.mat
%init_mat=zeros(size(con_mat,1)*size(triplet,1),size(pairs,1));
init_mat2=[];

for i=1:size(triplet,1)
%%
con_mat=(matcon_gen2ff(np_array(1,triplet(i,:))));
%con_mat=(matcon_gen2(np_array(1,triplet(i,:))));

if(i==1)
part_12=con_mat{1};
part_13=con_mat{2};
part_23=con_mat{3};
end

clear con_mat
% part_zeros=zeros(size(part_23));

%%




init_mat1=[];
ct_pairs=nchoosek(triplet(i,:),2);
ct_pairs_merged=str2num(strcat(num_str_balance(ct_pairs(:,1),length(np_array)),num_str_balance(ct_pairs(:,2),length(np_array))));
[loct,~]=find(pairs_merged==ct_pairs_merged');
count=1;
    for j=1:size(pairs,1)
        if(((sum(loct==j))~=0) && count==1)
            temp_mat=(part_12);
            count=count+1;
            %clear part_12
             
        elseif((sum(loct==j))~=0 && count==2)
            temp_mat=(part_13);
            count=count+1;
            %clear part_13

        
        
        elseif((sum(loct==j))~=0 && count==3)
            temp_mat=part_23;
             count=1;
             %clear part_23
        end
        if((sum(loct==j))==0 )
            part_zeros=sparse(prod(np_array(1,triplet(i,:))),prod(np_array(1,pairs(j,:))));
            temp_mat=(part_zeros);
            %clear part zeros;
            
            %count=1;
        end
        
       init_mat1=([init_mat1 temp_mat]); 
        
    end
init_mat2=([init_mat2;init_mat1]);    
    
end

A_ineqportion=(init_mat2);




%%



end