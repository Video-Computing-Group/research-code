function out= matcon_gen2ff(np_array)
m=np_array(1);
n=np_array(2);
p=np_array(3);


temp1=ones(p,1);
use_12=repmat(temp1,1,m*n);
[r12,c12] = size(use_12);
i12     = 1:numel(use_12);
j12   = repmat(1:c12,r12,1);
BLOCK_12     = sparse(i12',j12(:),use_12(:));



[row_1,col_1]=find(BLOCK_12);
coeff_vec=(0:m:m*p-m)';
col_1m=rem(col_1,m)+m*(max((1-rem(col_1,m)),0));
col_2=col_1m+repmat(coeff_vec,[m*n,1]);
% BLOCK_13=sparse(row_1,col_2,1,size(BLOCK_12,1),m*p);
BLOCK_13=sparse(row_1,col_2,1);

[~,col_12]=ind2sub([m n],col_1);
[~,col_13]=ind2sub([m p],col_2);
linearInd23_col=sub2ind([n p],col_12,col_13);
% BLOCK_23=sparse(row_1,linearInd23_col,1,size(BLOCK_12,1),n*p);
BLOCK_23=sparse(row_1,linearInd23_col,1);

out{1}=(BLOCK_12);
out{2}=(BLOCK_13);
out{3}=(BLOCK_23);

      
