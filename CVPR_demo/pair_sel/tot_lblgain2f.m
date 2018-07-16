function [tot_gain,lmat_struct]= tot_lblgain2f(tlmatstruct,lable_collect,cam_n)
camNamearr=1:cam_n;
triplet=nchoosek(camNamearr,3);
pairs=nchoosek(camNamearr,2); 
pairs_merged=str2num(strcat(num2str(pairs(:,1)),num2str(pairs(:,2))));
for i=1:size(pairs,1) 
ct_cama=pairs(i,1);
ct_camb=pairs(i,2);
[rowa, ~]=find(triplet==ct_cama); 
[rowb, ~]=find(triplet==ct_camb); 
common_rows=intersect(rowa,rowb); 
lmat_ct=lable_collect{i};
ct_campair=pairs_merged(i); 
for j=1:length(common_rows)
    %label_temp=zeros(size(lmat_ct)); %%% temp. lable holder 
    ct_row=triplet(common_rows(j),:); 
    ct_pairs=nchoosek(ct_row,2); 
    ct_pairs_merged=str2num(strcat(num2str(ct_pairs(:,1)),num2str(ct_pairs(:,2)))); 
    tran_pairs=ct_pairs_merged(find(ct_pairs_merged~=ct_campair));
    [row ,~]=find(pairs_merged==tran_pairs');  
    pair_1=pairs(row(1),:);
    pair_2=pairs(row(2),:);
    pair_1label=lable_collect{row(1)}; 
    pair_2label=lable_collect{row(2)}; 
   
    pair_1labelt=tlmatstruct{row(1)}; 
    pair_2labelt=tlmatstruct{row(2)}; 


    common_point=intersect(pair_1,pair_2); 
    if (pair_1(1)==common_point)
     pair_1label=pair_1label';
     pair_1labelt=pair_1labelt';
    end
        
    if (pair_2(1)~=common_point)
     pair_2label=pair_2label';
     pair_2labelt=pair_2labelt';
    end
%%%%%   
%This portion extracts labels using eq. relations 
tempx1=pair_1label.*pair_1labelt;
tempx2=pair_2label.*pair_2labelt;
lt1=logical(tempx1*pair_2label);
lt2=logical(pair_1label*tempx2);
label_temp=or(lt1,lt2);

% label_temp=logical(pair_1label*pair_2label);

% %label_notavl_temp=not(logical(double(not(pair_1labelt))*double(not(pair_2labelt))));
%  label_notavl_temp=(double(not(tempx1))*double(not(tempx2)));
% label_notavl_temp(label_notavl_temp<size(pair_1label,2))=0;
% label_temp=not(logical(label_notavl_temp));
% label_temp=logical(tempx1*tempx2);
%label_temp=lt3;
% sum(sum(label_temp))-sum(sum(sanity2))

% label_temp=label_temp.*label_available;
% sanity_check2=sum(sum(label_temp~=sanity2))
%%%%%
    
    
    
    
    lmat_ct=double(or(lmat_ct,label_temp));
end
lmat_struct{i}=lmat_ct;
ct_tlmat=tlmatstruct{i};
tot_gain(i)=sum(sum(ct_tlmat.*lmat_ct))/sum(sum(ct_tlmat));

    
end
% bias=bias=



end

function [out1,out2 ]= sanity_check(pair_1label,pair_2label,pair_1labelt,pair_2labelt)
pair_1label=pair_1label';
pair_1labelt=pair_1labelt';
out1=zeros(size(pair_1label,2),size(pair_2label,2));
out2=zeros(size(pair_1label,2),size(pair_2label,2));
for i=1:size(pair_1label,1)

     [col_tran1]=find(pair_1label(i,:));
    [col_tran2]=find(pair_2label(i,:));
    out1(col_tran1,col_tran2)=1;   
    
 for nec_1=1:length(col_tran1)
     for nec_2=1:length(col_tran2)
        if(pair_1labelt(i,col_tran1(nec_1))==1 ||pair_2labelt(i,col_tran2(nec_2))==1)
           out2(col_tran1(nec_1),col_tran2(nec_2))=1;
        end
     end
 end

    
end


end


