%% Voxel-wise comparison between PC-SPGR and PC-bSSFP
%
% - script generates histograms or AP/FH/LR/magnitude velocity values
% within the pipes of the flow phantom
% - FlowPhantom9 scans used for analysis - these were matched resolutions
% and FOVs
%
%


%% Data Choices
qflowDataDir = 'data_raw_rigid_to_ap';
bssfpDataDir = 'data_rigid_to_qflow';

numBssfpStacks = 5;
extVelDir = '_noSag';


%% Admin
fp9Dir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9';

qflowReconDir = [fp9Dir '\QFlow'];
bssfpReconDir = [fp9Dir '\PC-bSSFP'];


%% Load QFlow data
cd(qflowReconDir);
cd(qflowDataDir);

switch qflowDataDir
        
        % QFlow from .nii
        case 'data'

            velFiles = dir('*_vel.nii.gz');
            % velFiles = dir('*_vel_*scaled.nii.gz');
            velFiles = {velFiles([1,2,3]).name};
            for ii = 1:numel(velFiles)    
                nii = load_nii(velFiles{ii});   
                VEL(:,:,:,ii) = nii.img;
            end

        % QFlow from .nii, with FFD registration
        case 'data_ffd'
            
            velFiles = dir('*_vel_*scaled.nii.gz');
            velFiles = {velFiles([1,2,3]).name};
            for ii = 1:numel(velFiles)    
                nii = load_nii(velFiles{ii});   
                VEL(:,:,:,ii) = nii.img;
            end
            
        % QFlow from .raw, with no registration
        case 'data_raw'

            phFiles = dir('*_ph_corr.nii.gz');
            phFiles = {phFiles([1,2,3]).name};
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, with rigid registration
        case 'data_raw_rigid'
            
            phFiles = dir('*_ph_corr_rigid.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, with FFD registration
        case 'data_raw_ffd'

            phFiles = dir('*_ph_corr_ffd.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, rigid registered to QFlow AP stack
        case 'data_raw_rigid_to_ap'
            
            phFiles = dir('*_ph_corr_rigid.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);
            
        % QFlow from .raw, affine registered to QFlow AP stack
        case 'data_raw_affine_to_ap'
        
            phFiles = dir('*_ph_corr_affine.nii.gz');
            phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
            for ii = 1:numel(phFiles)    
                nii = load_nii(phFiles{ii});   
                PH(:,:,:,ii) = nii.img;
            end
            venc = 30; VEL = venc * (PH./pi);       
        
end

vq = double(VEL);



%% Load PC-bSSFP data
cd(bssfpReconDir);
cd(bssfpDataDir);
cd([num2str(numBssfpStacks) '_phase_stacks_recon' extVelDir]);

switch bssfpDataDir
        
        % PC-bSSFP transformed to axial stack
        case 'data_trans'
            
            nii = load_nii('ap_vel_trans.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_trans.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_trans.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
        % PC-bSSFP rigid registered to axial bSSFP stack
        case 'data_rigid'
            
            nii = load_nii('ap_vel_rigid.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_rigid.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_rigid.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
        % PC-bSSFP FFD registered to axial bSSFP stack
        case 'data_ffd'
            
            nii = load_nii('ap_vel_ffd.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_ffd.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_ffd.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
        % PC-bSSFP rigid registered to axial QFlow stack
        case 'data_rigid_to_qflow'
            
            nii = load_nii('ap_vel_rigid.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_rigid.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_rigid.nii.gz');
            vb(:,:,:,3) = double(nii.img);


        % PC-bSSFP affine registered to axial QFlow stack
        case 'data_affine_to_qflow'

            nii = load_nii('ap_vel_affine.nii.gz');
            vb(:,:,:,1) = double(nii.img);
            nii = load_nii('fh_vel_affine.nii.gz');
            vb(:,:,:,2) = -1.*double(nii.img);
            nii = load_nii('rl_vel_affine.nii.gz');
            vb(:,:,:,3) = double(nii.img);
            
            
end


%% Create Velocity Magnitude Volumes
vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% Load Water Pipe Mask
cd(qflowReconDir);
cd mask

P = load_nii('mask/inside_pipes_mask3.nii.gz');
% S = load_nii('spherical_part_mask.nii.gz');

M = double(P.img);


%% View Velocity Images
implay_RR([M.*vq(:,:,:,1), M.*vq(:,:,:,2), M.*vq(:,:,:,3), M.*vq3d;...
           M.*vb(:,:,:,1), M.*vb(:,:,:,2), M.*vb(:,:,:,3), M.*vb3d]...
           ,'jet',[-10,10]);     
       

%% Plot Histograms

lowerThresh = -20;
upperThresh = 20;
binSize = 0.2;

bkrdClip = 0.025;

for ii = 1:3
    M_vq(:,:,:,ii) = M.*vq(:,:,:,ii);
    M_vb(:,:,:,ii) = M.*vb(:,:,:,ii);
end

M_vq3d = M.*vq3d;
M_vb3d = M.*vb3d;

% Set background/really low values to NaN
% - NaNs excluded from histogram function
% - use QFlow velocity magnitude for thresholding
nanMask = double(M_vq3d>bkrdClip); nanMask(nanMask<1)=NaN;

for ii = 1:3
    tempVect = M_vq(:,:,:,ii).*nanMask; QVhist(:,ii) = tempVect(:);
    tempVect = M_vb(:,:,:,ii).*nanMask; BVhist(:,ii) = tempVect(:);
end

tempVect = M_vq3d .* nanMask; QVhist(:,4) = tempVect(:);
tempVect = M_vb3d .* nanMask; BVhist(:,4) = tempVect(:);


%% Components
figure;
for ii = 1:3
    
    subplot(2,2,ii);
    histogram(QVhist(:,ii),'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');
    hold on;
    histogram(BVhist(:,ii),'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');
    
    legend('PC-SPGR','PC-bSSFP');
    xlabel('Velocity [cm/s]');
    ylabel('Counts');
    

    if ii == 1
        title('AP');
    elseif ii == 2
        title('FH');
    elseif ii == 3
        title('RL');
    end
end


%% Magnitude
subplot(2,2,4);
histogram(QVhist(:,4),'binEdges',0:binSize/2:upperThresh,'Normalization','Probability');
hold on;
histogram(BVhist(:,4),'binEdges',0:binSize/2:upperThresh,'Normalization','Probability');

legend('PC-SPGR','PC-bSSFP');
xlabel('Velocity [cm/s]');
ylabel('Counts');
title('Magnitude');


%% Old stuff for PDFs
% figure; hold on;
% Qpdca = fitdist(QVhist(QVhist<10),'Normal');
% Bpdca = fitdist(BVhist(BVhist<10),'Normal');

% Qpdca = fitdist(QVhist,'Normal');
% Bpdca = fitdist(BVhist,'Normal');
% 
% figure; hold on;
% plot([0:0.25:10],pdf(Qpdca,[0:0.25:10]));
% plot([0:0.25:10],pdf(Bpdca,[0:0.25:10]));
% legend('PC-SPGR','PC-bSSFP');
% xlabel('Velocity [cm/s]');


%% For Bland-Altmans
BA.ap(:,1)  = QVhist(~isnan(BVhist(:,1))); %ap
BA.ap(:,2)  = BVhist(~isnan(BVhist(:,1)));
BA.fh(:,1)  = QVhist(~isnan(BVhist(:,2))); %fh
BA.fh(:,2)  = BVhist(~isnan(BVhist(:,2)));
BA.rl(:,1)  = QVhist(~isnan(BVhist(:,3))); %rl
BA.rl(:,2)  = BVhist(~isnan(BVhist(:,3)));
BA.mag(:,1) = QVhist(~isnan(BVhist(:,4))); %magnitude
BA.mag(:,2) = BVhist(~isnan(BVhist(:,4)));

