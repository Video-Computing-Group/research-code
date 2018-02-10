%% Code Author: Ahmed Tashrif Kamal - tashrifahmed@gmail.com
% http://www.ee.ucr.edu/~akamal
% no permission necessary for non-commercial use
% Date: 4/27/2013

%% 
disp('Generating Observations...');

[zt,zIdt,zCountt] = observeWithClutter(lamda,Nc,Nt,xa,H,R,FOV,T);

save ('savedObservations','zt','zIdt','zCountt');