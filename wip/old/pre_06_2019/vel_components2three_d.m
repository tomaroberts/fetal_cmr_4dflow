function vel_components2three_d(nii1,nii2,nii3,outputFileName)

% quick function to convert 3 velocity components into a single
% 4- or 5-dimensional velocity volume. Useful for MRtrix

if nargin < 4
    outputFileName = 'vel3D.nii.gz';
end

v1 = load_untouch_nii(nii1);
v2 = load_untouch_nii(nii2);
v3 = load_untouch_nii(nii3);

dataDims = numel(size(v1.img));

%% save 4D-velocity volume
v3D = v1; %initalise nii structure

if dataDims == 3
    v3D.img = cat(dataDims+1,v1.img(:,:,:,1), v2.img(:,:,:,1), v3.img(:,:,:,1));
end

if dataDims == 4
    v3D.img = cat(dataDims+1,v1.img(:,:,:,:,1), v2.img(:,:,:,:,1), v3.img(:,:,:,:,1));
end

v3D.hdr.dime.dim(dataDims+2) = 3;
v3D.hdr.dime.pixdim(dataDims+2) = 1;
save_untouch_nii(v3D,outputFileName);

disp(['Saved 4-D velocity volume as: ' outputFileName ' ...']);

end