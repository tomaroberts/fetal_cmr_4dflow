sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP_MIRTK_recon';
cd(sp);

dataDir = '\data';
ktReconDir = '\ktrecon';
maskDir = '\mask';
velDir = '\vel_vol';
magDir = '\mag_vol';


%% load
cv = load_untouch_nii([sp magDir '\mag_vol.nii.gz']);
v0 = load_untouch_nii([sp velDir '\velocity-0.nii.gz']);
v1 = load_untouch_nii([sp velDir '\velocity-1.nii.gz']);
v2 = load_untouch_nii([sp velDir '\velocity-2.nii.gz']);

% v4D.img = sqrt( v0.img.^2 + v1.img.^2 + v2.img.^2 );
% 
% implay_RR([cv.img(:,:,:,1)]);  % slices
% implay_RR([cv.img(:,:,18,:)]); % time
% 
% implay_RR([v0.img(:,:,:,1), v1.img(:,:,:,1), v2.img(:,:,:,1)],'jet');     % slices
% implay_RR([v0.img(:,:,18,:), v1.img(:,:,18,:), v2.img(:,:,18,:)],'jet');  % time
% 
% implay_RR(v4D.img(:,:,18,:));  % time


%% save 3D-velocity volumes for MRTrix

v3 = v0;
v3.img = cat(4,v0.img(:,:,:,1), v1.img(:,:,:,1), v2.img(:,:,:,1));
v3.hdr.dime.dim(1) = 4;
v3.hdr.dime.dim(5) = 3;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_v0v1v2.nii.gz']);

v3 = v0;
v3.img = cat(4,v0.img(:,:,:,1), v2.img(:,:,:,1), v1.img(:,:,:,1));
v3.hdr.dime.dim(1) = 4;
v3.hdr.dime.dim(5) = 3;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_v0v2v1.nii.gz']);

v3 = v0;
v3.img = cat(4,v1.img(:,:,:,1), v0.img(:,:,:,1), v2.img(:,:,:,1));
v3.hdr.dime.dim(1) = 4;
v3.hdr.dime.dim(5) = 3;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_v1v0v2.nii.gz']);

v3 = v0;
v3.img = cat(4,v1.img(:,:,:,1), v2.img(:,:,:,1), v0.img(:,:,:,1));
v3.hdr.dime.dim(1) = 4;
v3.hdr.dime.dim(5) = 3;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_v1v2v0.nii.gz']);

v3 = v0;
v3.img = cat(4,v2.img(:,:,:,1), v0.img(:,:,:,1), v1.img(:,:,:,1));
v3.hdr.dime.dim(1) = 4;
v3.hdr.dime.dim(5) = 3;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_v2v0v1.nii.gz']);

v3 = v0;
v3.img = cat(4,v2.img(:,:,:,1), v1.img(:,:,:,1), v0.img(:,:,:,1));
v3.hdr.dime.dim(1) = 4;
v3.hdr.dime.dim(5) = 3;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_v2v1v0.nii.gz']);