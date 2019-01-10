%% Time Starts
tic

%% Read and Normalize the image
tempVar  = load('./reuse_variables.mat');
dir = tempVar.directory;
extension = '.tif';
imageInd = tempVar.imageInd;
imageToRead = strcat(dir,imageInd{3});
im = readImage(imageToRead);
% thresholdVec =[0.0775 0.0550 0.0750 0.0500 0.0650 0.0575 0.0650 0.0650 0.0750 0.0400 0.0550 0.0800];
thresholdVec = thresholdVec;%   0.0625    0.0550    0.0400    0.0200];
% thresholdVec =[0.0450   0.0500  0.095   0.075   0.075   0.075   0.065   0.095   0.09    0.060   0.050];
sizeIm = size(im);
Masks = zeros(sizeIm(1),sizeIm(2),length(imageInd)-2); 
clear BGImages;
for k=1:length(imageInd)-2
    filename = strcat(dir,imageInd{k+2});
    tempIm = readImage(filename);
    BGImages(:,:,k) = tempIm;
    gBlurIm = imfilter(tempIm,fspecial('gaussian',[10 10],15));
    [B cellStructureTemp MASK1] = subSegmentation(imhmin(gBlurIm,thresholdVec(k)),k);
    numCells{k} = length(B);
    Masks(:,:,k) = MASK1;
    cellStructure{k} = cellStructureTemp;
%     ShowOverlayMask(tempIm,MASK1, k);
    phi(:,:,k) = getPhi(tempIm,MASK1);  
%     k
end
save 'BGImages.mat' BGImages;
save 'numCells.mat' numCells;
save 'Masks.mat' Masks;
save 'cellStructure.mat' cellStructure;
save 'phi.mat' phi;

%% Time Stops
toc
