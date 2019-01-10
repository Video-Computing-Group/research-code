% HIT F5 to run this file
%% Time Starts
clear;clc;close all;
tic
delete *.mat
%%Path to folder containing images
directory = './example_input/';
file_to_save_variables = './reuse_variables.mat';
%% Read and Normalize the image
%extension = '.tif';
%imageInd = [1 5];
listDir = dir(directory);
imageInd = extractfield(listDir,'name');
threshRange = [0.02 0.08];
thresholdJump = 0.005;

readIm = strcat(directory,imageInd{3});
im = readImage(readIm);
sizeIm = size(im);

threshMetrics = zeros(length(imageInd),...
    int16((threshRange(2)-threshRange(1))/thresholdJump) + 1,4); 

for k=1:length(imageInd)-2
    filename = strcat(directory,imageInd{k+2});
    tempIm = readImage(filename);
    gBlurIm = imfilter(tempIm,fspecial('gaussian',[10 10],15));
    count=1;
    for threshold = threshRange(1):thresholdJump:threshRange(2)
        MASK1 = watershed(imhmin(gBlurIm,threshold));
        MASK1 = MASK1==0;
        MASK1 = bwareaopen(MASK1,100);
        MASK1 = double(MASK1);
        [B L] = bwboundaries(MASK1, 'holes');       
        [centroid area cellPixelArray] = findCellInfo(B,L);
        threshMetrics(k,count,:)=[threshold length(area), mean(area), std(area)];
        count=count+1;
    end
end
save('threshStudyMetrics', 'threshMetrics');
save(file_to_save_variables,'directory');
save(file_to_save_variables,'imageInd','-append');

%% Time Stops
toc
threshStudy;
getSegment;
%figure,imshow(phi(:,:,1))