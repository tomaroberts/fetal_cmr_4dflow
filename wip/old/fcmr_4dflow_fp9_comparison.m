%% Voxel-wise comparison between PC-SPGR and PC-bSSFP
%
% - script generates histograms or AP/FH/LR/magnitude velocity values
% within the pipes of the flow phantom
% - FlowPhantom9 scans used for analysis - these were matched resolutions
% and FOVs
%
% - REQUIRES: 2019_07_04_Flow_Phantom9_PC-bSSFP_*.m to generate ph/velocity
% files.
%
%


%% Data Choices
% qflowDataDir = 'data_raw_rigid_to_ap';
% bssfpDataDir = 'data_rigid_to_qflow';

qflowDataDir = 'data_raw_affine_to_ap';
bssfpDataDir = 'data_affine_to_qflow';

numBssfpStacks = 5;
extVelDir = '_noSag';

maskFileName = 'inside_pipes_mask3.nii.gz';


%% Admin
fp9Dir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9';

qflowReconDir = [fp9Dir '\QFlow'];
bssfpReconDir = [fp9Dir '\PC-bSSFP'];


%% Load QFlow data
cd(qflowReconDir);
cd(qflowDataDir);

MAGq = load_nii('ap_mag_affine.nii.gz'); MAGq = MAGq.img;

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

MAGb = load_nii('ax_mag_affine.nii.gz'); MAGb = MAGb.img;

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

P = load_nii(['mask/' maskFileName]);
% S = load_nii('spherical_part_mask.nii.gz');

M = double(P.img);


%% View Velocity Images
% p = [1,2,3];
% d = 2; %180 degrees
% 
% 
% implay_RR([rot_ms(permute(M.*vq(:,:,:,1),p),d), rot_ms(permute(M.*vq(:,:,:,2),p),d), rot_ms(permute(M.*vq(:,:,:,3),p),d), rot_ms(permute(M.*vq3d(:,:,:,1),p),d);...
%            rot_ms(permute(M.*vb(:,:,:,1),p),d), rot_ms(permute(M.*vb(:,:,:,2),p),d), rot_ms(permute(M.*vb(:,:,:,3),p),d), rot_ms(permute(M.*vb3d(:,:,:,1),p),d)]...
%            ,'jet',[-10,10]); 
%        
% implay_RR([rot_ms(permute(MAGq,p),d); rot_ms(permute(MAGb,p),d)],'gray'); 
     

%% Images for paper: 
% Mold = M;
% M = M(5:165,9:169,:);
implay_RR([rot_ms(permute(M.*vq(5:165,9:169,:,1),p),d), rot_ms(permute(M.*vq(5:165,9:169,:,2),p),d), rot_ms(permute(M.*vq(5:165,9:169,:,3),p),d), rot_ms(permute(M.*vq3d(5:165,9:169,:,1),p),d);...
           rot_ms(permute(M.*vb(5:165,9:169,:,1),p),d), rot_ms(permute(M.*vb(5:165,9:169,:,2),p),d), rot_ms(permute(M.*vb(5:165,9:169,:,3),p),d), rot_ms(permute(M.*vb3d(5:165,9:169,:,1),p),d)]...
           ,'jet',[-10,10]);

d = 1; %180 degrees
rcrop = 11:108;
ccrop = 7:104;
sl = 61;

Mold = M;
M = M(rcrop,ccrop,sl);

imtar([fliplr(rot90(MAGq(rcrop,ccrop,sl),d)); fliplr(rot90(MAGb(rcrop,ccrop,sl),d))]);
xticklabels(''); yticklabels(''); xticks(''); yticks(''); colorbar('off');

% imtar([rot90(M.*vq(rcrop,ccrop,sl,1),d), rot90(M.*vq(rcrop,ccrop,sl,2),d), rot90(M.*vq(rcrop,ccrop,sl,3),d), rot90(M.*vq3d(rcrop,ccrop,sl,1),d);...
%        rot90(M.*vb(rcrop,ccrop,sl,1),d), rot90(M.*vb(rcrop,ccrop,sl,2),d), rot90(M.*vb(rcrop,ccrop,sl,3),d), rot90(M.*vb3d(rcrop,ccrop,sl,1),d)]...
%        ,-10,10); colormap('jet'); xticklabels(''); yticklabels(''); xticks(''); yticks(''); colorbar('off'); set(gca,'Visible','off')

cd('C:\Users\tr17\Documents\Papers\4D_Velocity_CINE_Paper\MRM_v1');
imtar(fliplr(rot90(M.*vq(rcrop,ccrop,sl,1),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-spgr_Vap.tiff');
imtar(fliplr(rot90(M.*vq(rcrop,ccrop,sl,2),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-spgr_Vfh.tiff');
imtar(fliplr(rot90(M.*vq(rcrop,ccrop,sl,3),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-spgr_Vrl.tiff');
imtar(fliplr(rot90(M.*vq3d(rcrop,ccrop,sl,1),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-spgr_Vmagnitude.tiff');

imtar(fliplr(rot90(M.*vb(rcrop,ccrop,sl,1),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-bssfp_Vap.tiff');
imtar(fliplr(rot90(M.*vb(rcrop,ccrop,sl,2),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-bssfp_Vfh.tiff');
imtar(fliplr(rot90(M.*vb(rcrop,ccrop,sl,3),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-bssfp_Vrl.tiff');
imtar(fliplr(rot90(M.*vb3d(rcrop,ccrop,sl,1),d)),-10,10,'paperFig'); colormap('jet'); saveas(gcf,'FP9_pc-bssfp_Vmagnitude.tiff');

imtar(fliplr(rot90(MAGq(rcrop,ccrop,sl),d)),[],[],'paperFig'); colormap('gray'); saveas(gcf,'FP9_pc-spgr_Structural.tiff');


%% Plot Histograms

lowerThresh = -20;
upperThresh = 20;
binSize = 0.25;

bkrdClip = 0;

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
% - plotted/biases calculated in GraphPad file:
% E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\FP9_BlandAltman_affine_to_qflow.pzfx
clear BA

BA.ap(:,1)  = rmmissing(QVhist(:,1));
BA.ap(:,2)  = rmmissing(BVhist(:,1));
BA.ap(:,3)  = mean(BA.ap,2);
BA.ap(:,4)  = 100*( (BA.ap(:,1)-BA.ap(:,2)) ./ BA.ap(:,3) );
% BA.ap(find(abs(BA.ap(:,4)) > spuriousVal),:) = [];
BA.ap(find(sum(sign(BA.ap(:,1:2)),2) == 0 ),:) = [];

BA.fh(:,1)  = rmmissing(QVhist(:,2));
BA.fh(:,2)  = rmmissing(BVhist(:,2));
BA.fh(:,3)  = mean(BA.fh,2);
BA.fh(:,4)  = 100*( (BA.fh(:,1)-BA.fh(:,2)) ./ BA.fh(:,3) );
% BA.fh(find(abs(BA.fh(:,4)) > spuriousVal),:) = [];
BA.fh(find(sum(sign(BA.fh(:,1:2)),2) == 0 ),:) = [];

BA.rl(:,1)  = rmmissing(QVhist(:,3));
BA.rl(:,2)  = rmmissing(BVhist(:,3));
BA.rl(:,3)  = mean(BA.rl,2);
BA.rl(:,4)  = 100*( (BA.rl(:,1)-BA.rl(:,2)) ./ BA.rl(:,3) );
% BA.rl(find(abs(BA.rl(:,4)) > spuriousVal),:) = [];
BA.rl(find(sum(sign(BA.rl(:,1:2)),2) == 0 ),:) = [];

BA.mag(:,1) = rmmissing(QVhist(:,4));
BA.mag(:,2) = rmmissing(BVhist(:,4));
BA.mag(:,3)  = mean(BA.mag,2);
BA.mag(:,4)  = 100*( (BA.mag(:,1)-BA.mag(:,2)) ./ BA.mag(:,3) );





