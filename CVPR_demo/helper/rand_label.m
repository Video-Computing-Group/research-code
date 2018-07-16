function a= rand_label(r,c,n)
a=zeros(r,c);
index = randperm(numel(a));
a(index(1:n)) = 1;
a=a';



end