function out =greedy_selff2(f,A,lim)
[~,ind]=sort(f);
ind_greedy=zeros(size(f));
A=A(:,ind);

%A_t=A';
i=1;
while sum(ind_greedy)<lim

if(i<3)
 ind_greedy(i)=1;   
    
else
 ind_greedy(i)=1;
 flag=find(A*ind_greedy>2);
 
 if(sum(flag)>0)
   ind_greedy(i)=0;  
 end
 
 
end
if(i==14)
end

 i=i+1;
 if(i>length(f))
     break
 end
 

 
end

out=zeros(size(f));
out(ind(find(ind_greedy)))=1;
    
    
end



