% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Truncates each member functions
function F = trunc(F,c)
for i = 1:length(F.Fs)
	F.Fs{i}=sfo_fn_trunc(F.Fs{i},c);
end
