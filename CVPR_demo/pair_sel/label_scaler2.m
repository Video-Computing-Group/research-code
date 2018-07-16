function B=label_scaler2(A)
B=A;
B(A==0)=-1;
B(A>0)=1;



end