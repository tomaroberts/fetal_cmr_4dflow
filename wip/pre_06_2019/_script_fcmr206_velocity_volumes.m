sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr206';
cd(sp);

dataDir = '\data';
ktReconDir = '\ktrecon';
maskDir = '\mask';
velDir = '\vel_vol';
cineDir = '\cine_vol';


%% load
cv = load_untouch_nii([sp cineDir '\cine_vol.nii.gz']);
v0 = load_untouch_nii([sp velDir '\velocity-0.nii.gz']);
v1 = load_untouch_nii([sp velDir '\velocity-1.nii.gz']);
v2 = load_untouch_nii([sp velDir '\velocity-2.nii.gz']);

v4D.img = sqrt( v0.img.^2 + v1.img.^2 + v2.img.^2 );

implay_RR([cv.img(:,:,:,1)]);  % slices
implay_RR([cv.img(:,:,18,:)]); % time

implay_RR([v0.img(:,:,:,1), v1.img(:,:,:,1), v2.img(:,:,:,1)],'jet');     % slices
implay_RR([v0.img(:,:,18,:), v1.img(:,:,18,:), v2.img(:,:,18,:)],'jet');  % time

implay_RR(v4D.img(:,:,18,:));  % time


%% save 3D-velocity volume for MRTrix
v3 = v0;
v3.img = cat(4,v0.img(:,:,:,1), v1.img(:,:,:,1), v2.img(:,:,:,1));
v3.hdr.dime.dim(5) = 3;
v3.hdr.dime.pixdim(5) = 1;
save_untouch_nii(v3,[sp velDir '\vel3D_mrtrix_vol_t1.nii.gz']);