function out= matcon_gen2(np_array)
m=np_array(1);
n=np_array(2);
p=np_array(3);


temp1=ones(p,1);
blkp=flip(dec2bin(m*n));
blk_ct=[];
temp2=[];
for i=1:length(blkp)
    if(i==1)
    blk_ct=blkdiag(temp1,blk_ct);
    else
    blk_ct=blkdiag(blk_ct,blk_ct);
    end
    
    
   if(str2num(blkp(i))==1)
        blk_f=blkdiag(temp2,blk_ct);
        temp2=blk_f;
        
        
    end
    
      
    
end
BLOCK_12=sparse(blk_f);

BLOCK_13=sparse(size(BLOCK_12,1),m*p);


[row_1,col_1]=find(BLOCK_12);
coeff_vec=(0:m:m*p-m)';
col_1m=rem(col_1,m)+m*(max((1-rem(col_1,m)),0));
col_2=col_1m+repmat(coeff_vec,[m*n,1]);
linearInd_13 = sub2ind(size(BLOCK_13),row_1,col_2);
BLOCK_13(linearInd_13')=1;


%[~,col_2]=find(BLOCK_13);

BLOCK_23=sparse(size(BLOCK_12,1),n*p);
[~,col_12]=ind2sub([m n],col_1);
[~,col_13]=ind2sub([m p],col_2);
linearInd23_col=sub2ind([n p],col_12,col_13);
linearInd_23 =sub2ind(size(BLOCK_23),row_1,linearInd23_col);
BLOCK_23(linearInd_23')=1;

out{1}=sparse(BLOCK_12);
out{2}=sparse(BLOCK_13);
out{3}=sparse(BLOCK_23);

      
