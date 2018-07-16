function lb_gain_ratio = lbg_ratio(Lmat_opt,master_budget)

%%% computes how much label gain we have achieved 
tot=0;
tot_size=0;
for i=1:length(Lmat_opt)
    A=full(Lmat_opt{1,i});
    A(A==0)=1;
    temp=length(find(A==1));
    tot=tot+temp;
    tot_size=tot_size+length(A(:));
    
    
    
end


%lb_gain_ratio= ((tot-master_budget)/tot_size);
lb_gain_ratio= ((tot)/tot_size);



end