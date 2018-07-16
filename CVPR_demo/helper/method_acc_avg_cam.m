function [out] = method_acc_avg_cam(Method)
%%% calculates accuracy over all cameras
length(Method.full_set.acc);
num_pairs=length(Method.full_set.acc);
num_methods=4;
Method.methods{6}='Full Set';
Method.methods{3}='Half-Approx';


for i=1:num_methods
    temp=[];
for j=1:num_pairs
   temp=[temp;Method.acc{1,i}{j,1}(1)];   
end
out(i)=100*mean(temp);
end
temp=[];
for j=1:num_pairs
   temp=[temp;Method.full_set.acc{j,1}(1)];   
end
out(num_methods+1)=100*mean(temp);




% legend(Method.methods{1},Method.methods{2},Method.methods{3},Method.methods{4},Method.methods{6})
end
