% QFlow from .nii
% 2mm isotropic voxels
%
% Had to resample the 2mm .nii to match the bSSFP data. Did this with mirtk


%% Step 1)
%% separate magnitude and magnitude2 images in ISDPACS nifti files

% reconDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_04_Flow_Phantom8\QFlow';
% cd(reconDir);
% cd data_2mm
% 
% nii = load_untouch_nii('*FH*it0.nii.gz');
% nii.img = nii.img(:,:,1:2:end);
% nii.hdr.dime.dim(4) = size(nii.img,3);
% save_untouch_nii(nii,'fh_mag.nii.gz');
% 
% nii = load_untouch_nii('*RL*it0.nii.gz');
% nii.img = nii.img(:,:,1:2:end);
% nii.hdr.dime.dim(4) = size(nii.img,3);
% save_untouch_nii(nii,'rl_mag.nii.gz');
% 
% nii = load_untouch_nii('*AP*it0.nii.gz');
% nii.img = nii.img(:,:,1:2:end);
% nii.hdr.dime.dim(4) = size(nii.img,3);
% save_untouch_nii(nii,'ap_mag.nii.gz');

% *it3.nii.gz files renamed as fh_vel.nii.gz, etc.

%% Step 2)
%% Resample to bSSFP resolution
%e.g: mirtk resample-image fh_vel.nii.gz fh_vel_resample.nii.gz -imsize 176 176 120 -size 0.8523 0.8523 1.25 -interp bspline


%% Step 3)
%% Load stacks and scale
reconDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_04_Flow_Phantom8\QFlow';
cd(reconDir);
cd data_2mm

magFiles = dir('*_mag.nii.gz');
magFiles = {magFiles([1,2,3]).name}; %re-order: ap/fh/rl

velFiles = dir('*_vel_resample_unscaled.nii.gz');
velFiles = {velFiles([1,2,3]).name};

for ii = 1:numel(magFiles)
    nii = load_untouch_nii(magFiles{ii});
    MAG(:,:,:,ii) = double(nii.img);
    
    nii = load_untouch_nii(velFiles{ii});
    VEL(:,:,:,ii) = double(nii.img);
    
    % taken from nii header:
    slope = 0.0147;
    inter = -30;
    VEL(:,:,:,ii) = ( VEL(:,:,:,ii) .* slope ) + inter;
    
%     nii.img = VEL(:,:,:,ii);
%     nii.hdr.dime.datatype = 16;
%     save_untouch_nii(nii,[velFiles{ii}(1:end-16) '.nii.gz']);
    
end

clear nii

permOrder = [2,1,3];

% implay_RR([permute(MAG(:,:,:,1),permOrder), ...
%            permute(MAG(:,:,:,2),permOrder), ...
%            permute(MAG(:,:,:,3),permOrder)] ...
%            ,'gray');

implay_RR([permute(VEL(:,:,:,1),permOrder), ...
           permute(VEL(:,:,:,2),permOrder), ...
           permute(VEL(:,:,:,3),permOrder)] ...
           ,'jet',[-5,5]);

       