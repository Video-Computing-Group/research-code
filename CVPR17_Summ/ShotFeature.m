%% Code to find shot level feature mapping from C3D features
% Input: C3D feature file generated from ConcateC3D.m code
% Output: Final feature mat file
%--------------------------------------------------------------------------
clc; clear all;close all;
path = 'Y:\Rameswar\ECCV_2016\Features\CoSum\C3DFeatures';
path1 = 'Y:\Rameswar\Datasets\Single-view\CoSum\ShotSegmentation';
root_feature_dir_info = dir(path);
for k = 4:12
    feature_path = strcat(path,'\',root_feature_dir_info(k).name); 
    frm_num_path = strcat(path1,'\',root_feature_dir_info(k).name);
    feature_dir_info = dir([feature_path '\*.mat']);
    frm_num_dir_info = dir([frm_num_path '\*.mat']);
    feature_fileNames = {feature_dir_info.name};
    frm_num_fileNames = {frm_num_dir_info.name};
    nfeature = length (feature_fileNames);
    nfrm_num = length(frm_num_fileNames);
    for j = 1:nfrm_num
        feature_folder_path = strcat(feature_path,'\',feature_dir_info(j).name); 
        frm_num_folder_path = strcat(frm_num_path,'\',frm_num_dir_info(j).name); 
        clear C F frm_num feature_mat
        load(feature_folder_path);
        load(frm_num_folder_path);
        C = transpose(frm_num);
        [~,Nf]=size(F);
        track_prev=1;
        feature_mat=zeros(4096,length(C)/2);
        flag=1;
        for i=1:2:length(C)
            t = C(i+1)/16; u = floor(t);
            if (u+1) > Nf
                u1 = Nf-1;
            else
                u1 = u;
            end
            u = u1;
            if track_prev < Nf
            if (t - u) >= 0.5 % will always satify -- always true
                    v =  sum(transpose(F(:,track_prev:u+1)));
                    w =((u+1) - track_prev) + 1;
                    mean_v = v./w;
                    track_prev=u+2;
                else
                    v = sum(transpose(F(:,track_prev:u)));
                    w = (u - track_prev) + 1;
                    mean_v = v./w;
                    track_prev=u+1;
            end
            feature_mat(:,flag)=mean_v;
            flag=flag+1;        
            end
        end
        feature_mat=feature_mat(:,any(feature_mat)); % Removing zero columns
        s = strcat('C3D','_',feature_dir_info(j).name);
        save([feature_path,'\',s],'feature_mat');
    end
end