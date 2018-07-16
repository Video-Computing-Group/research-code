function [tot_gain,lmat_struct]= tot_lblgain(tlmatstruct,lable_collect,cam_n)
camNamearr=1:cam_n;
triplet=nchoosek(camNamearr,3);
pairs=nchoosek(camNamearr,2); %% [1-2; 1-3; 1-4; 2-3; 2-4; 3-4]
pairs_merged=str2num(strcat(num2str(pairs(:,1)),num2str(pairs(:,2))));
for i=1:size(pairs,1)
ct_cama=pairs(i,1);
ct_camb=pairs(i,2);
[rowa, ~]=find(triplet==ct_cama);
[rowb, ~]=find(triplet==ct_camb);   
common_rows=intersect(rowa,rowb);
lmat_ct=lable_collect{i};
ct_campair=pairs_merged(i);
% num_pcam1a=p_countpairwise(i,1);
% num_pcam1b=p_countpairwise(i,2);
% lmat_init=reshape(selvec(bias+1:bias+num_pcam1a*num_pcam1b,1),num_pcam1a,num_pcam1b);

for j=1:length(common_rows)
    label_temp=zeros(size(lmat_ct));
    ct_row=triplet(common_rows(j),:);
    ct_pairs=nchoosek(ct_row,2);
    ct_pairs_merged=str2num(strcat(num2str(ct_pairs(:,1)),num2str(ct_pairs(:,2))));
    tran_pairs=ct_pairs_merged(find(ct_pairs_merged~=ct_campair));
    [row ,~]=find(pairs_merged==tran_pairs');  %%%%%% check 
    pair_1=pairs(row(1),:);
    pair_2=pairs(row(2),:);
    pair_1label=lable_collect{row(1)};
    pair_2label=lable_collect{row(2)};


    common_point=intersect(pair_1,pair_2);
    if (pair_1(1)~=common_point)
     pair_1label=pair_1label';
    end
        
    if (pair_2(1)~=common_point)
     pair_2label=pair_2label';
    end
    
%     [~,col_tran1]=find(pair_1label);
%     [~,col_tran2]=find(pair_2label);
%     
%     label_temp(col_tran1,col_tran2)=1;
for loop_n=1:size(label_temp,1)
    [col_tran1]=find(pair_1label(loop_n,:));
    [col_tran2]=find(pair_2label(loop_n,:));
    label_temp([col_tran1],[col_tran])=1;

    
    
    
    
end
   
    
    
    
    
    
    
    
    
    lmat_ct=double(or(lmat_ct,label_temp));
end
lmat_struct{i}=lmat_ct;
ct_tlmat=tlmatstruct{i};
tot_gain(i)=sum(sum(ct_tlmat.*lmat_ct))/sum(sum(ct_tlmat));

    
end
% bias=bias=



end