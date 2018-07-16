function [tot_gain,lmat_struct,overhead,L_collect]= tot_lblgain3f(L_collect,tlmatstruct,lable_collect,simstruct,cam_n,reg,st_need,reg_main)

if nargin < 8
    reg_main=8;
end


if nargin < 7
    reg_main=8;
    st_need=20;
end

camNamearr=1:cam_n;
triplet=nchoosek(camNamearr,3);
pairs=nchoosek(camNamearr,2); %% [1-2; 1-3; 1-4; 2-3; 2-4; 3-4]
pairs_merged=str2num(strcat(num2str(pairs(:,1)),num2str(pairs(:,2))));
overhead=0;
reg_num=1+1/(reg_main*reg); %% 8*reg
if (st_need <15 && reg>350)
 reg_num=1-1/(4*reg); %% 8*reg   
end

if reg>449
    reg_num=-10;
    
end
%  if(if reg_num<
for i=1:size(pairs,1)
ct_cama=pairs(i,1);
ct_camb=pairs(i,2);
[rowa, ~]=find(triplet==ct_cama);
[rowb, ~]=find(triplet==ct_camb);   
common_rows=intersect(rowa,rowb); %%% find which triplets of cameras contain both current camera a and b 
lmat_ct=lable_collect{i}; %%  available labels of current pair (a,b)
sim_ct=simstruct{i};  %% access the similarity matrix for current pair
ct_campair=pairs_merged(i);
% num_pcam1a=p_countpairwise(i,1);
% num_pcam1b=p_countpairwise(i,2);
% lmat_init=reshape(selvec(bias+1:bias+num_pcam1a*num_pcam1b,1),num_pcam1a,num_pcam1b);

for j=1:length(common_rows) %% iterate over relevant triplets for (a,b)
    label_temp=zeros(size(lmat_ct));
    ct_row=triplet(common_rows(j),:); %% current triplet
    ct_pairs=nchoosek(ct_row,2);
    ct_pairs_merged=str2num(strcat(num2str(ct_pairs(:,1)),num2str(ct_pairs(:,2))));
    tran_pairs=ct_pairs_merged(find(ct_pairs_merged~=ct_campair));
    [row ,~]=find(pairs_merged==tran_pairs');  %%%%%% check 
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
    
%This portion extracts labels using eq. relations 
tempx1=pair_1label.*pair_1labelt;
tempx2=pair_2label.*pair_2labelt;
lt1=logical(tempx1*pair_2label);
lt2=logical(pair_1label*tempx2);
label_temp=or(lt1,lt2);
lmat_ct=double(or(lmat_ct,label_temp));    
    
end
sim_ct_temp=sim_ct.*not(lmat_ct);
sim_ct_tempxx=sim_ct_temp>reg_num*(mean(sim_ct,2));
reg_sim=full(sim_ct_tempxx);
if(sum(reg_sim(:))>st_need)
 sim_spread=sim_ct_temp.*sim_ct_tempxx;
 sim_spread=sim_spread(:);
 spread_l=zeros(size(sim_spread));
 [~,ind_spread]=sort(sim_spread,'descend');
 ind_st=ind_spread(1:st_need);
 spread_l(ind_st)=1;
 sim_ct_tempxx=reshape(spread_l,size(sim_ct_tempxx));
end
% lmat_ct=double(or(lmat_ct,sim_ct_temp));
lmat_ct=double(or(lmat_ct,sim_ct_tempxx));
L_collect{i}=sparse(double(or(L_collect{i},sim_ct_tempxx)));

overhead=overhead+sum(sum(sim_ct_tempxx));
st_need=st_need-sum(sum(sim_ct_tempxx));
lmat_struct{i}=lmat_ct;
ct_tlmat=tlmatstruct{i};
tot_gain(i)=sum(sum(ct_tlmat.*lmat_ct))/sum(sum(ct_tlmat));

    
end
% bias=bias=

overhead;

end