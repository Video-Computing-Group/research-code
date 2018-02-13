% Helper function that makes an sfo_fn_* object Octave ready
% This is required since Octave handles function objects differently
%
% Author: Andreas Krause (krausea@gmail.com)
%
% function F = sfo_octavize(F_input)
% F_input: sfo_fn_* function object
%
% Example: F = sfo_octavize(sfo_fn_example);
%
function F = sfo_octavize(F_input)
F = sfo_fn_wrapper(@(A) F_input(A));
