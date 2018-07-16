function out = num_str_balance(a,bl_num)
max_digit_num=numel(num2str(bl_num));


for i=1:length(a)
    
    ct_digit_num=numel(num2str(a(i)));
    gap=max_digit_num-ct_digit_num;
    if(gap>0)
    out(i,:)= sprintf(sprintf('%%0%is',max_digit_num),num2str(a(i)));
    else
    out(i,:)=num2str(a(i));
    end
    
    
    
    
end





end