function [im] = readImage(imgPath)
im = imread(imgPath);
if length(size(im)) == 3
    try 
        if size(im,3) == 3
            im = rgb2gray(im);
        end
    catch
        im = im;
    end
end
%im = imresize(im,[450 450]);
im = normalize(double(im));

