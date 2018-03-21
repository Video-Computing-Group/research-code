%% Code to concate the C3D features extracted using Caffe
% It requires the read_binary_blob function given in C3D feature extraction folder
% Input: Folder containing all the fc6-1 files
% Output: Feature Mat file 
%------------------------------------------------------------------------------------
clc;clear all;close all;
Path = 'Y:\Rameswar\ECCV_2016\Features\TVSUM50\C3DFeatures\BT\'; % Change path accordingly
dir_info = dir(Path);
for i = 3:7
FolderPath = strcat(Path,dir_info(i).name); 
filelist = dir([FolderPath '\*.fc6-1']);
fileNames = {filelist.name};
nFeatures = length(fileNames);
F = zeros(4096,nFeatures);
for j=1:nFeatures
    [s,data] = read_binary_blob([FolderPath '/' fileNames{j}]);
    F(:,j) = data';
end
save(FolderPath,'F');
end

