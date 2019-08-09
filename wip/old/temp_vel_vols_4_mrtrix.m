%% Fetal data
%- ie: temporal component

v0 = load_untouch_nii('velocity-final-0.nii.gz');
v1 = load_untouch_nii('velocity-final-1.nii.gz');
v2 = load_untouch_nii('velocity-final-2.nii.gz');


% DEFAULT (ie: no flipping)   

% saveDir = 'vel_vol_4d/default/';
% saveName = '_velocity-final-4D-t';


% SINGLE FLIP    

% saveDir = 'vel_vol_4d/0flipped/';
% saveName = '_0flipped-velocity-final-4D-t';
% v0.img = -1 .* v0.img;

% saveDir = 'vel_vol_4d/1flipped/';
% saveName = '_1flipped-velocity-final-4D-t';
% v1.img = -1 .* v1.img;

saveDir = 'vel_vol_4d/2flipped/'; %%% WORKS IN FLOW PHANTOM
saveName = '_2flipped-velocity-final-4D-t';
v2.img = -1 .* v2.img;

% TWO FLIPS

% saveDir = 'vel_vol_4d/0flipped1flipped/';
% saveName = '_0flipped1flipped-velocity-final-4D-t';
% v0.img = -1 .* v0.img;
% v1.img = -1 .* v1.img;

% saveDir = 'vel_vol_4d/0flipped2flipped/';
% saveName = '_0flipped2flipped-velocity-final-4D-t';
% v0.img = -1 .* v0.img;
% v2.img = -1 .* v2.img;

% saveDir = 'vel_vol_4d/1flipped2flipped/'; %%% POSSIBLY LOOKS BETTER THAN JUST 2flipped
% saveName = '_1flipped2flipped-velocity-final-4D-t';
% v1.img = -1 .* v1.img;
% v2.img = -1 .* v2.img;

% ALL FLIPS

% saveDir = 'vel_vol_4d/0flipped1flipped2flipped/';
% saveName = '_0flipped1flipped2flipped-velocity-final-4D-t';
% v0.img = -1 .* v0.img;
% v1.img = -1 .* v1.img;
% v2.img = -1 .* v2.img;


% saveDir = '';
mkdir(saveDir);



%% save 5D volume

% 
% v3D = v0;
% v3D.img = cat(5,v0.img, v1.img, v2.img);
% v3D.img = permute(v3D.img,[1,2,3,5,4]);
% 
% v3D.hdr.dime.dim(5) = size(v3D.img,4);
% v3D.hdr.dime.dim(6) = size(v3D.img,5);
% 
% v3D.hdr.dime.pixdim(6) = v3D.hdr.dime.pixdim(5);
% v3D.hdr.dime.pixdim(5) = 1;
% 
% save_untouch_nii(v3D,[saveDir 'v3D_2flipped.nii.gz']);



for ii = 1:size(v0.img,4)
% for ii = 20
    
    %v012
    v3D = v0;
    v3D.img = cat(4,v0.img(:,:,:,ii), v1.img(:,:,:,ii), v2.img(:,:,:,ii));
    
    v3D.hdr.dime.dim(5) = 3;
    v3D.hdr.dime.pixdim(5) = 1;
    save_untouch_nii(v3D,[saveDir 'v012' saveName num2str(ii) '.nii.gz']);
    
%     %v021
%     v3D = v0;
%     v3D.img = cat(4,v0.img(:,:,:,ii), v2.img(:,:,:,ii), v1.img(:,:,:,ii));
%     
%     v3D.hdr.dime.dim(5) = 3;
%     v3D.hdr.dime.pixdim(5) = 1;
%     save_untouch_nii(v3D,[saveDir 'v021' saveName num2str(ii) '.nii.gz']);
    
    
    
    
    
%     %v102
%     v3D = v0;
%     v3D.img = cat(4,v1.img(:,:,:,ii), v0.img(:,:,:,ii), v2.img(:,:,:,ii));
%     
%     v3D.hdr.dime.dim(5) = 3;
%     v3D.hdr.dime.pixdim(5) = 1;
%     save_untouch_nii(v3D,[saveDir 'v102' saveName num2str(ii) '.nii.gz']);
%     
%     %v120
%     v3D = v0;
%     v3D.img = cat(4,v1.img(:,:,:,ii), v2.img(:,:,:,ii), v0.img(:,:,:,ii));
%     
%     v3D.hdr.dime.dim(5) = 3;
%     v3D.hdr.dime.pixdim(5) = 1;
%     save_untouch_nii(v3D,[saveDir 'v120' saveName num2str(ii) '.nii.gz']);
% 
%     
%     
%     %v201
%     v3D = v0;
%     v3D.img = cat(4,v2.img(:,:,:,ii), v0.img(:,:,:,ii), v1.img(:,:,:,ii));
%     
%     v3D.hdr.dime.dim(5) = 3;
%     v3D.hdr.dime.pixdim(5) = 1;
%     save_untouch_nii(v3D,[saveDir 'v201' saveName num2str(ii) '.nii.gz']);
% 
%     %v210
%     v3D = v0;
%     v3D.img = cat(4,v2.img(:,:,:,ii), v1.img(:,:,:,ii), v0.img(:,:,:,ii));
%     
%     v3D.hdr.dime.dim(5) = 3;
%     v3D.hdr.dime.pixdim(5) = 1;
%     save_untouch_nii(v3D,[saveDir 'v210' saveName num2str(ii) '.nii.gz']);

    
end

disp('DONE ...');









%% Phantom data
%- ie: single timepoint

vv=load_untouch_nii('velocity-vector-final.nii.gz');
% vv=load_untouch_nii('velocity-vector-final_cm_per_sec.nii.gz');

% saveDir = 'vel_vol_4d/default/';

saveDir = 'vel_vol_4d/flipped3/';
vv.img(:,:,:,3) = vv.img(:,:,:,3) .* -1;

mkdir(saveDir);

v123=vv;

vv.img = cat(4,v123.img(:,:,:,1), v123.img(:,:,:,2), v123.img(:,:,:,3));
save_untouch_nii(vv,[saveDir 'velocity-vector-final-123.nii.gz']);

vv.img = cat(4,v123.img(:,:,:,1), v123.img(:,:,:,3), v123.img(:,:,:,2));
save_untouch_nii(vv,[saveDir 'velocity-vector-final-132.nii.gz']);

vv.img = cat(4,v123.img(:,:,:,2), v123.img(:,:,:,1), v123.img(:,:,:,3));
save_untouch_nii(vv,[saveDir 'velocity-vector-final-213.nii.gz']);

vv.img = cat(4,v123.img(:,:,:,2), v123.img(:,:,:,3), v123.img(:,:,:,1));
save_untouch_nii(vv,[saveDir 'velocity-vector-final-231.nii.gz']);

vv.img = cat(4,v123.img(:,:,:,3), v123.img(:,:,:,1), v123.img(:,:,:,2));
save_untouch_nii(vv,[saveDir 'velocity-vector-final-312.nii.gz']);

vv.img = cat(4,v123.img(:,:,:,3), v123.img(:,:,:,2), v123.img(:,:,:,1));
save_untouch_nii(vv,[saveDir 'velocity-vector-final-321.nii.gz']);



