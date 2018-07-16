function out =greedy_sel(f,A,lim)
[~,ind]=sort(f);
ind_greedy=zeros(size(f));
A=A(:,ind);
i=1;
while sum(ind_greedy)<lim

if(i<3)
 ind_greedy(i)=1;   
    
else
   
 temp1=find(A(:,i)); %%% find current rows involved   
 temp2=(A(temp1,:)); %% find other columns involved for those rows
 temp2(:,i)=[]; %%% current common column removal
 [rowt,colt]=find(temp2); %%% find other non-zero columns, they are column sorted 
 [~,rowind]=sort(rowt); %%% we have to row-sort them in r1-r1 , r2-r2, r3-r3 format r1<r2<r3 
 col_s=colt(rowind); %%% get the cols correspoding to the sorted rows
 col_s=col_s+1;
 temp3=ind_greedy(col_s);
 temp4=temp3(1:2:end).*temp3(2:2:end);
  if(sum(temp4)==0)
  ind_greedy(i)=1;       
  end
 flag=find(A*ind_greedy>2);
 if(i==13)
     pause
 end
 
 if(sum(flag)>0)
 end
  
end
 
 i=i+1;
 if(i>length(f))
     break
 end
end

out=zeros(size(f));
out(ind(find(ind_greedy)))=1;
    
    
end



