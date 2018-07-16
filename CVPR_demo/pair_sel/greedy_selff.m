function out =greedy_selff(f,A,lim)
[~,ind]=sort(f);
ind_greedy=zeros(size(f));
A=A(:,ind);
tic
[row_main,col_main]=find(A);
A_t=sparse(col_main,row_main,1);
toc
%A_t=A';
i=1;
while sum(ind_greedy)<lim

if(i<3)
 ind_greedy(i)=1;   
    
else
 ind_greedy(i)=1;       
  
 [temp1]=find(A(:,i)); 
 [rx,~]=find(A_t(:,temp1));

temp3=ind_greedy(rx);
temp4=temp3(1:3:end).*temp3(2:3:end).*temp3(3:3:end);
  if(sum(temp4)~=0)
  ind_greedy(i)=0;       
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



