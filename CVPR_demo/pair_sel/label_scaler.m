
function scaled_d=label_scaler(A)
scaled_d=A;
scaled_d(A>0)=1;
scaled_d(A<0)=-1;
%scaled_d=A;


end





