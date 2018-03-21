%% Code for Temporal Segmentation of a Video (Dividing videos into shots)
% Reference Paper: Video Co-Sum, CVPR'15
% Idea: Measure the amount of changes (sum-squared pixel-wise difference) 
% between 2 consecutive frames in both RGB and HSV color spaces.
%% Load the video
clc; clear all; close all
Path = 'Y:\Rameswar\Datasets\Single-view\TVSUM50\frames\VU'; % Input is a folder containing videos of the same topic
dir_info = dir(Path);
for j = 3:7
    FolderPath = strcat(Path,'\',dir_info(j).name); 
    filelist = dir([FolderPath '/*.jpg']);
    fileNames = {filelist.name};
    nImages = length(fileNames);

    diff_rgb = zeros(nImages-1,1);
    diff_hsv = zeros(nImages-1,1);
    for k = 1:1:(nImages-1)
        input_I1 = imread([FolderPath '/' fileNames{k}]);
        input_I2 = imread([FolderPath '/' fileNames{k+1}]);
        I1_rgb = im2double(input_I1);
        I2_rgb = im2double(input_I2);
        I1_vec_rgb = I1_rgb(:);
        I2_vec_rgb = I2_rgb(:);
        I1_hsv = rgb2hsv(input_I1);
        I2_hsv = rgb2hsv(input_I2);
        I1_vec_hsv = I1_hsv(:);
        I2_vec_hsv = I2_hsv(:);
        diff_rgb(k) = sum((I1_vec_rgb - I2_vec_rgb).^2);
        diff_hsv(k) = sum((I1_vec_hsv - I2_vec_hsv).^2);  
    end
    total_diff = diff_rgb + diff_hsv;
    flag = 2;
    frm_num(1)= 1;
    temp = frm_num(1);
    for l = 1:(nImages-2)
        if (((abs(total_diff(l) - total_diff(l+1))/total_diff(l))) > 0.75)
            if (abs(temp - (l+1)) >= 32 && abs(temp - (l+1)) <= 96)
                frm_num(flag)= l+1;
                frm_num(flag+1) = l+1;
                flag = flag+2;
                temp = l+1;
            else
                while (abs(temp - (l+1)) > 96 )
                    frm_num(flag) = temp + 96;
                    frm_num(flag+1) = temp + 96;
                    flag = flag + 2;
                    temp = temp + 96;
                end
            end
        end
    end
        frm_num(flag)= nImages;
        s = strcat('frm_num','_',dir_info(j).name);
        save(s,'frm_num');
        clear frm_num
end