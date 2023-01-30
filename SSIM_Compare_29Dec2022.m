%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSIM_Compare.m
% Written by Jinxia (Fiona) Yao and Bradley Fitzgerald
% Last modified: 29 December 2022

%%%% Script Description %%%%
% This script takes two brain delay maps (intended for comparing CO2-MRI
% with DSC-MRI) and compares them using the structure element of the 
% structural similarity index measure (SSIM; see Zhou Wang 2014)

%%%% Primary Inputs %%%%%%%%%%%%%%%%%%

% DSC_img = blood arrival delay map computed via DSC-MRI

% CO2_img = blood arrival delay map computed via CO2-MRI

% GM_mask = gray matter mask

% WM_mask = white matter mask


%%%% Primary Outputs %%%%%%%%%%%%%%%%%

% structure_map = map of voxel-wise structural values from SSIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

%% Input images to compare (just for testing, not in final version)
% First input CO2 and DSC delay maps
% We assume that all empty voxels have NaN value
addpath('/work/BF/Yunjie/Functions/NIfTI_20140122/');
addpath('/work/BF/Yunjie/Functions');
DSC_nifti=load_nii_gz( '/work/BF/Yunjie/Brads_Versions/CO2_Paper/Scripts_2022/mean_DSC_standard_3_shift.nii.gz');
CO2_nifti=load_nii_gz( '/work/BF/Yunjie/Brads_Versions/CO2_Paper/Scripts_2022/ALL_SUBS_Averaged_CO2_Delays_20Dec2021.nii.gz');

DSC_img = DSC_nifti.img;
CO2_img = CO2_nifti.img;

% Input GM and WM masks
filename = '/work/BF/templates/MNIparcs/MNI152_T1_2mm_brain_seg.nii.gz';
segmap_nii = load_nii_gz(filename);
GM_mask = segmap_nii.img==2;
WM_mask = segmap_nii.img==3;

%% Input CO2 and DSC delay maps
% Replace these values with your input CO2 and DSC delay maps
DSC_img = DSC_img;
CO2_img = CO2_img;

%% Input gray matter (GM) and white matter (WM) masks
% Replace these values with your input GM and WM masks
GM_mask = GM_mask;
WM_mask = WM_mask;

%% Create image mask which includes valid voxels from both images
dsc_img_mask = ~isnan(DSC_img);
filled_dsc_img_mask = logical(imfill(double(dsc_img_mask)));
empty_voxels = dsc_img_mask ~= filled_dsc_img_mask;
DSC_empty_mask_indx = find(empty_voxels);

img_mask_co2 = ~isnan(CO2_img);
final_mask = dsc_img_mask .*img_mask_co2;

%% Demean DSC and CO2 delay maps
DSC_map = DSC_img - nanmean(reshape(DSC_img, [], 1));
mean_dsc = nanmean(reshape(DSC_map, [], 1));
DSC_map(~final_mask) = mean_dsc;

CO2_map = CO2_img;
mean_co2 = nanmean(reshape(CO2_map, [], 1));
CO2_map(~final_mask) = mean_co2;

%% Compute SSIM map based only on structure element
[score1, structure_map] = ssim(DSC_map, CO2_map, 'Exponents', [0 0 1]);
structure_map(DSC_empty_mask_indx) = nan; 
structure_map(isnan(CO2_img)) = nan;
structure_map(isnan(DSC_img)) = nan;

%% Compute metrics

%Average structure element over whole brain
str_score = mean(reshape(structure_map,[],1),'omitnan');
disp(strcat('Average structure element over whole brain: ', num2str(str_score)));

%Percentage of voxels with high structural similarity
data = reshape(structure_map, [], 1);
data(isnan(data)) = [];
perc_high = nnz(data>0.5) / nnz(data);
disp(strcat('Percentage of voxels with high structural similarity:', num2str(perc_high)));

%% Compute GM and WM structural similarity metrics
str_GM = structure_map;
str_GM(GM_mask==0) = NaN;
str_WM = structure_map;
str_WM(WM_mask==0) = NaN;

mean_GM = nanmean(reshape(str_GM, [], 1))
mean_WM = nanmean(reshape(str_WM, [], 1))

perc_high_GM = nnz(str_GM>0.5) / nnz(~isnan(str_GM))
perc_high_WM = nnz(str_WM>0.5) / nnz(~isnan(str_WM))

