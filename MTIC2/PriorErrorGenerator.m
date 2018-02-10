%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
disp('Generating Prior Errors...');
eta = cell(Nt,1);
for j=1:Nt
    eta{j} = mvnrnd(zeros(p,1),P)';
end
save('PriorError.mat','eta');
