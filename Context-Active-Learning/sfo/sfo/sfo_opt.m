% Andreas Krause (krausea@cs.cmu.edu)
% Generate a default option struct used to specify options in various
% algorithms.
%
% function opt = sfo_opt(vals)
% vals: cell array of form vals = {'NAME1',VAL1,'NAME2',VAL2,...}
%
% Example: opt = sfo_opt({'val',3.14,'color','red'})

function opt = sfo_opt(vals)
opt = []; 
opt.binsearch_tolerance = 1e-5;
opt.verbosity_level = 1;
if exist('vals','var')
    for i = 1:2:length(vals)
        opt = setfield(opt,vals{i},vals{i+1});
    end
end
