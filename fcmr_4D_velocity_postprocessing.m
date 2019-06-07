%% 4D Cardiac Flow - fcmr post-processing
%
% - Applies polynomial correction to velocity volumes to clean up 'static'
% areas
% - Writes .vtk files for usage in Paraview
% - Writes .nii files for usage in MRtrix
% --- Both correctly formatted so velocity vectors point in correct
% direction
%
%


%% Admin
fNum = 194;
% fNum = 202;
% fNum = 254;
fcmrDir = ['C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr' num2str(fNum)];
cd(fcmrDir);

velDir = 'vel_vol';


%% Run velocity volume polynomial correction
% fcmr_pc_velocity_correction( fcmrDir , velDir );


%% Load cine volume
cd('cine_vol')
if ~isfile('cine_vol-RESLICE.nii.gz')
    reslice_nii('cine_vol.nii.gz', 'cine_vol-RESLICE.nii.gz');
end
cine_nii = load_nii('cine_vol-RESLICE.nii.gz');

nFrame = size(cine_nii.img,4);


%% Load velocity volume
cd(['../' velDir]);

% %%% straight out of SVRTK
% if ~isfile('velocity-final-RESLICE-0.nii.gz') || ~isfile('velocity-final-RESLICE-1.nii.gz') || ~isfile('velocity-final-RESLICE-2.nii.gz') 
%     reslice_nii('velocity-final-0.nii.gz', 'velocity-final-RESLICE-0.nii.gz');
%     reslice_nii('velocity-final-1.nii.gz', 'velocity-final-RESLICE-1.nii.gz');
%     reslice_nii('velocity-final-2.nii.gz', 'velocity-final-RESLICE-2.nii.gz');
% end
% 
% velx_nii = load_nii('velocity-final-RESLICE-0.nii.gz');
% vely_nii = load_nii('velocity-final-RESLICE-1.nii.gz');
% velz_nii = load_nii('velocity-final-RESLICE-2.nii.gz');


%%% velocity volume polyCorr version
if ~isfile('velocity-final-polyCorr-RESLICE-0.nii.gz') || ~isfile('velocity-final-polyCorr-RESLICE-1.nii.gz') || ~isfile('velocity-final-polyCorr-RESLICE-2.nii.gz') 
    reslice_nii('velocity-final-polyCorr-0.nii.gz', 'velocity-final-polyCorr-RESLICE-0.nii.gz');
    reslice_nii('velocity-final-polyCorr-1.nii.gz', 'velocity-final-polyCorr-RESLICE-1.nii.gz');
    reslice_nii('velocity-final-polyCorr-2.nii.gz', 'velocity-final-polyCorr-RESLICE-2.nii.gz');
end

velx_nii = load_nii('velocity-final-polyCorr-RESLICE-0.nii.gz');
vely_nii = load_nii('velocity-final-polyCorr-RESLICE-1.nii.gz');
velz_nii = load_nii('velocity-final-polyCorr-RESLICE-2.nii.gz');


%% Load masks

cd('../mask');

% Admin - blood pool mask:
% BP mask automatically created in fcmr_pc_velocity_correction.m
maskBloodPoolFileName = 'mask_blood_pool';

% Admin - custom masks:
if fNum == 194
    % All:
%     masksToUse = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC'};
    % Right-sided
%     masksToUse = {'mask_RV','mask_ROT','mask_LA_RA','mask_IVC_SVC'};
    % Left-sided
%     masksToUse = {'mask_aorta','mask_LV','mask_LOT'};
    % Left-sided w/out LOT
    masksToUse = {'mask_aorta','mask_LV'};
    %     maskBloodPoolFileName = 'mask_blood_pool_cleaned';
elseif fNum == 202
    masksToUse = {'mask_aorta','mask_LV','mask_RV','mask_LOT','mask_ROT','mask_LA_RA','mask_IVC_SVC','mask_PA_DA'};
end



% Load blood pool mask
% open mask / reslice if necessary
if ~isfile([maskBloodPoolFileName '-RESLICE.nii.gz'])
    reslice_nii([maskBloodPoolFileName '.nii.gz'], [maskBloodPoolFileName '-RESLICE.nii.gz']);
end

mask_bp = load_nii([maskBloodPoolFileName '-RESLICE.nii.gz']);
mask_bp.img = double(mask_bp.img);
for ii = 1:nFrame; mask_bp.img(:,:,:,ii) = mask_bp.img(:,:,:,1); end


% Load custom masks
for mm = 1:numel(masksToUse)      
    
    % current mask
    maskFileName = masksToUse{mm};
    
    % open mask / reslice if necessary
    if ~isfile([maskFileName '-RESLICE.nii.gz'])
        reslice_nii([maskFileName '.nii.gz'], [maskFileName '-RESLICE.nii.gz']);
    end

    mask = load_nii([maskFileName '-RESLICE.nii.gz']);
    mask.img = double(mask.img);
    
    % fix mask_aorta/mask_IVC_SVC for fcmr194 (to do with extra voxels when re-slicing):
    if fNum == 194 && strcmp(maskFileName,'mask_aorta') == 1 || fNum == 194 && strcmp(maskFileName,'mask_IVC_SVC') == 1
        for tt = 1:nFrame; mask_re.img(:,:,:,tt) = imresize3(mask.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt))); end
        mask.img = mask_re.img; clear mask_re;
    end
    
    % initalise maskCombined
    if mm == 1
        maskCombined.img = zeros(size(mask.img));
    end

    maskCombined.img = maskCombined.img + mask.img;
%     maskCombined.img(:,:,1:64,:) = maskCombined.img(:,:,1:64,:) + mask.img(:,:,1:64,:);
    
end

mask.img = single(logical(maskCombined.img));

for ii = 1:nFrame; mask.img(:,:,:,ii) = mask.img(:,:,:,1); end


%% Make component velocity volumes for VTK
Vx = velx_nii.img;
Vy = vely_nii.img;
Vz = velz_nii.img;

% cm/s
Vx = 1e2 .* Vx;
Vy = 1e2 .* Vy; 
Vz = 1e2 .* Vz;


%% Resize
% FIXME: for some reason mask is 1 voxel larger than original volume...? 
for tt = 1:nFrame
    Vx_re(:,:,:,tt) = imresize3(Vx(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vy_re(:,:,:,tt) = imresize3(Vy(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    Vz_re(:,:,:,tt) = imresize3(Vz(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    mask_re.img(:,:,:,tt) = imresize3(mask.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
    mask_bp_re.img(:,:,:,tt) = imresize3(mask_bp.img(:,:,:,tt),size(cine_nii.img(:,:,:,tt)));
end

Vx = Vx_re; clear Vx_re;
Vy = Vy_re; clear Vy_re;
Vz = Vz_re; clear Vz_re;
mask.img = mask_re.img; clear mask_re;
mask_bp.img = mask_bp_re.img; clear mask_bp_re;


%% Apply blood pool mask to cine_nii / custom masks to velocity volumes

% magnitude
cine_masked_nii.img = double(cine_nii.img) .* mask_bp.img;

% velocity
Vx_masked = Vx .* mask.img;
Vy_masked = Vy .* mask.img;
Vz_masked = Vz .* mask.img;

Vmag_masked = sqrt(Vx_masked.^2 + Vy_masked.^2 + Vz_masked.^2);

% cm/s
% nb: QFlow already in cm/s


%% Make meshgrids for VTK
[x,y,z] = meshgrid(1:size(Vx_masked,2),1:size(Vx_masked,1),1:size(Vx_masked,3));


%% Eliminate spurious values from Vmag
% improves Paraview visualisation
% hist(nonzeros(Vmag_masked(:)))
% errIdx = find(Vmag_masked(:) > 50);
% Vx_masked(errIdx) = 0; Vy_masked(errIdx) = 0; Vz_masked(errIdx) = 0; Vmag_masked(errIdx) = 0;


%% Set Low Velocity Threshold
% - cleans up volumes for viewing / fewer vectors to render
lowerBoundIdx = find(Vmag_masked(:) < 1);
Vx_masked(lowerBoundIdx) = 0; Vy_masked(lowerBoundIdx) = 0; Vz_masked(lowerBoundIdx) = 0; Vmag_masked(lowerBoundIdx) = 0;


%% Write .vtk files
cd(fcmrDir);
cd(velDir);
% mkdir paraview
% cd paraview
% mkdir paraview_noPolyCorr
% cd paraview_noPolyCorr
mkdir paraview_left_side_noLOT
cd paraview_left_side_noLOT

for tt = 1:nFrame
    
    % magnitude
    vtkwrite(['cine_vol_masked_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'scalars', 'magnitude_intensity', cine_masked_nii.img(:,:,:,tt) );
         
%     %%% CORRECT CONFIGURATION OF VxVyVz in FLOW PHANTOM 6 ...
%     vtkwrite(['vel_vol_masked_VyVx-Vz_t' num2str(tt-1) '.vtk'], ...
%              'structured_grid', x, y, z, ... 
%              'vectors', 'vector_field', Vy_masked(:,:,:,tt), Vx_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt), ...
%              'scalars', 'velocity_magnitude', Vmag_masked(:,:,:,tt) );  
         
    %%% CORRECT CONFIGURATION FOR FETAL (advised by David)
    vtkwrite(['vel_vol_masked_VxVy-Vz_t' num2str(tt-1) '.vtk'], ...
             'structured_grid', y, x, z, ... 
             'vectors', 'vector_field', Vx_masked(:,:,:,tt), Vy_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt), ...
             'scalars', 'velocity_magnitude', Vmag_masked(:,:,:,tt) );  

end

disp('Finished making .vtk files ...');


%% Write .nii

cd(fcmrDir);
cd(velDir);
% mkdir mrtrix
% cd mrtrix
% mkdir mrtrix_noPolyCorr
% cd mrtrix_noPolyCorr
mkdir paraview_left_side_noLOT
cd paraview_left_side_noLOT

% get .nii to use as basis
Vx3D_nii = load_untouch_nii([fcmrDir '\' velDir '\velocity-final-polyCorr-RESLICE-0.nii.gz']);

for tt = 1:nFrame
    v3D = Vx3D_nii;
    v3D.img = cat(4,Vx_masked(:,:,:,tt), Vy_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt));
    
    v3D.hdr.dime.dim(1) = 4;
    v3D.hdr.dime.dim(2:5) = [size(Vx_masked,1) size(Vx_masked,2) size(Vx_masked,3) 3];
    v3D.hdr.dime.pixdim(5) = 0;
    
    save_untouch_nii(v3D,['VxVy-Vz_t' num2str(tt-1) '.nii.gz']);
end

disp('Finished making .nii files ...');


%% Write .vtk files with different lower bounds

% lbRange = [5,10,20,30];
% 
% cd(fcmrDir);
% cd(velDir);
% mkdir('paraview_lowerBounds');
% cd('paraview_lowerBounds');
% 
% for rr = 1:numel(lbRange)
%     
%     lowerBoundIdx = find(Vmag_masked(:) < lbRange(rr));
%     Vx_masked(lowerBoundIdx) = 0; Vy_masked(lowerBoundIdx) = 0; Vz_masked(lowerBoundIdx) = 0; Vmag_masked(lowerBoundIdx) = 0;
%     
%     for tt = 1:nFrame
%         
%         %%% CORRECT CONFIGURATION FOR FETAL (advised by David)
%         vtkwrite(['vel_lb' num2str(lbRange(rr)) ' t' num2str(tt-1) '.vtk'], ...
%              'structured_grid', y, x, z, ... 
%              'vectors', 'vector_field', Vx_masked(:,:,:,tt), Vy_masked(:,:,:,tt), -1 .* Vz_masked(:,:,:,tt), ...
%              'scalars', 'velocity_magnitude', Vmag_masked(:,:,:,tt) );  
%         
%         % rename lb5 files to lb05 for ease of use in Paraview:
%         if lbRange(rr) == 5
%             movefile( ['vel_lb' num2str(lbRange(rr)) ' t' num2str(tt-1) '.vtk'] , ['vel_lb0' num2str(lbRange(rr)) ' t' num2str(tt-1) '.vtk']);
%         end
%         
%     end
%     
%     disp(['Finished making .vtk files with lower bound = ' num2str(lbRange(rr)) ' ...']);
%     
% end


