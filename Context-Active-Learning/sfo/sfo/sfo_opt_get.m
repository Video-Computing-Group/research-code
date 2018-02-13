% Andreas Krause (krausea@cs.cmu.edu)
% Get a value from an option struct
%
% function val = sfo_opt_get(opt,name,default)
% opt: option struct (created using sfo_opt)
% name: field name
% default (optional): default value
%
% Example: TOL = sfo_opt_get(opt,'tolerance',1e-10)

function val = sfo_opt_get(opt,name,default)
if isfield(opt,name)
    val = getfield(opt,name);
else
    val = default;
end
