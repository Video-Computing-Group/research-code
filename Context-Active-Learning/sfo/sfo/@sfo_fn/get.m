% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function val = get(F, propName)
switch propName
case 'current_set'
   val = F.current_set;
case 'current_val'
   val = F.current_val;
otherwise
   error([propName,' Is not a valid asset property'])
end
