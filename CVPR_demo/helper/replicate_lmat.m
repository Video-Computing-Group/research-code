function rep_lmat = replicate_lmat(Lmat,Id1,Id2)
%%% replicate edges for multi-shot setting, Lmat is the single-shot label matrix,
%%% Id1 is the multi-shot ID vector in sorted identity order , Id2 is the
%%% multi-shot cam2 ID vector in sorted identity order 
%%% output: rep_lmat => replicated label matrix from single-shot label
%%% matrix lmat
cnt=1;
Id1_adjusted(1)=cnt;
for i=2:length(Id1)

  if (Id1(i)==Id1(i-1))
    Id1_adjusted(i)=cnt;
  else
    cnt=cnt+1;
    Id1_adjusted(i)=cnt;
      
end


end

cnt=1;
Id2_adjusted(1)=cnt;
for i=2:length(Id2)

  if (Id2(i)==Id2(i-1))
    Id2_adjusted(i)=cnt;
  else
    cnt=cnt+1;
    Id2_adjusted(i)=cnt;
      
end


end

A_temp=Lmat(Id1_adjusted,:);
rep_lmat=A_temp(:,Id2_adjusted);




end