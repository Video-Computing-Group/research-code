% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function a = set(a,varargin)
propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
   case 'current_val'
      a.current_val = val;
   case 'current_set'
      a.current_set = val;
   otherwise
      error('sfo_fn properties: current_val, current_set')
   end
end