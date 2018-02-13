% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information
%% sums utility functions across buckets
function val = inc(F,As,el)
As = sfo_unique_fast(As);
F = init(F,As);
m = length(F.Fs);
val = 0;
oldScore = get(F,'current_val');

if sum(As == el)>0
    val = oldScore;
    return
end

A_part = partition(F,As);
bucket = mod(el,m)+1;
sensor = floor(el/m);
for i = 1:m
    if i == bucket
        val = val + inc(F.Fs{i},A_part{i},sensor);
    else
        [F.Fs{i},v] = init(F.Fs{i},A_part{i});
        val = val + v;
%        val = val + get(F.Fs{i},'current_val');
    end
end
%F = set(F,'current_set',sfo_unique_fast([As el]),'current_val',val);
