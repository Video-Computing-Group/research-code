% Energy function for ising model for image denoising
% Implementation by Andreas Krause (krausea@gmail.com)
%
% function F = sfo_fn_ising(img,coeffPix,coeffH,coeffV,coeffD)
% img: n x m binary array (image)
% A: subset of the pixels set to 1 (ranging in 1: (n*m))
% coeffPix/H/V/Diag: negative log potentials for differing pixels,
%   horizontal, vertical or diagonal intensity mismatch
%
% Example: F = sfo_fn_ising(ones(3),1,1,1,0); F([1 4 7])

function F = sfo_fn_ising(img,coeffPix,coeffH,coeffV,coeffD)
val0 = evalIsing(img,[],coeffPix,coeffH,coeffV,coeffD);
fn = @(A) evalIsing(img,A,coeffPix,coeffH,coeffV,coeffD)-val0;
F = sfo_fn_wrapper(fn);

function E = evalIsing(img,A,coeffPix,coeffH,coeffV,coeffD)
A = sfo_unique_fast(A);
[r,c] = size(img);
mask = zeros(size(img));
mask(A) = 1;
delta = abs(mask-img);
Epix = sum(delta(:)); %energy through pixelwise disagreement
Ehor = abs(mask(:,2:c,:)-mask(:,1:(c-1))); Ehor = sum(Ehor(:)); %energy through "horizontal" pixel differences
Evert = abs(mask(2:r,:)-mask(1:(r-1),:)); Evert = sum(Evert(:)); %energy through "vertical" pixel differences
Ediag = abs(mask(1:(r-1),1:(c-1))-mask(2:r,2:c)) + abs(mask(2:(r),1:(c-1))-mask(1:(r-1),2:c)); Ediag = sum(Ediag(:));
E = (coeffPix*Epix+coeffH*Ehor(:)+coeffV*Evert(:)+coeffD*Ediag(:));
