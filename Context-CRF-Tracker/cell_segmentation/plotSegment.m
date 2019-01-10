function plotSegment(imNum)
tempVar  = load('reuse_variables');
directory = tempVar.directory;
tempVar2 = load('phi');
phi = tempVar2.phi;
imageInd = tempVar.imageInd;
imageToRead = strcat(directory,imageInd{imNum+2})
im = readImage(imageToRead);
%figure
try
    use_imNum = size(phi,3);
    %imshowpair(im,phi(:,:,imNum))
    figure, imshow(im)
    figure, imshow(phi(:,:,imNum))
catch
    imshowpair(im,phi(:,:))
    figure, imshow(im)
    figure, imshow(phi(:,:))
end
end



    
