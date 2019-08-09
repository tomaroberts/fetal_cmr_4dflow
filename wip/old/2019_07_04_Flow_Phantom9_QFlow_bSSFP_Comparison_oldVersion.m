%% qflow from nii (no FFD)
% qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data';
% cd(qDir);
% 
% velFiles = dir('*_vel.nii.gz');
% % velFiles = dir('*_vel_*scaled.nii.gz');
% velFiles = {velFiles([1,2,3]).name}; %re-order: ap/fh/rl
% 
% for ii = 1:numel(velFiles)    
%     nii = load_nii(velFiles{ii});   
%     VEL(:,:,:,ii) = nii.img;
% end
% 
% vq = double(VEL);
% 
% % mask
% M = load_nii('../mask/inside_pipes_mask.nii.gz');
% M = single(M.img);
% 
% % implay_RR([M.*vq(:,:,:,1), ...
% %            M.*vq(:,:,:,2), ...
% %            M.*vq(:,:,:,3)] ...
% %            ,'jet',[-5,5]);


%% qflow from nii with FFD
% qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_ffd';
% cd(qDir);
% 
% velFiles = dir('*_vel_*scaled.nii.gz');
% velFiles = {velFiles([1,2,3]).name}; %re-order: ap/fh/rl
% 
% for ii = 1:numel(velFiles)    
%     nii = load_nii(velFiles{ii});   
%     VEL(:,:,:,ii) = nii.img;
% end
% 
% vq = double(VEL);
% 
% % mask
% M = load_nii('../mask/inside_pipes_mask.nii.gz');
% M = single(M.img);
% 
% implay_RR([M.*vq(:,:,:,1), ...
%            M.*vq(:,:,:,2), ...
%            M.*vq(:,:,:,3)] ...
%            ,'jet',[-5,5]);


%% qflow from raw (no registration)
qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_raw';
cd(qDir);

phFiles = dir('*_ph_corr.nii.gz');
phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl

for ii = 1:numel(phFiles)    
    nii = load_nii(phFiles{ii});   
    PH(:,:,:,ii) = nii.img;
end

venc = 30;
VEL = venc * (PH./pi);

vq = double(VEL);

% mask
M = load_nii('../mask/inside_pipes_mask.nii.gz');
M = single(M.img);

% implay_RR([M.*vq(:,:,:,1), ...
%            M.*vq(:,:,:,2), ...
%            M.*vq(:,:,:,3)] ...
%            ,'jet',[-5,5]);


%% qflow from raw (rigid registered to qflow_ap)
qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_raw_rigid_to_ap';
cd(qDir);

phFiles = dir('*_ph_corr_rigid.nii.gz');
phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl

for ii = 1:numel(phFiles)    
    nii = load_nii(phFiles{ii});   
    PH(:,:,:,ii) = nii.img;
end

venc = 30;
VEL = venc * (PH./pi);

vq = double(VEL);

% mask
M = load_nii('../mask/inside_pipes_mask3.nii.gz');
M = single(M.img);

% implay_RR([M.*vq(:,:,:,1), ...
%            M.*vq(:,:,:,2), ...
%            M.*vq(:,:,:,3)] ...
%            ,'jet',[-5,5]);


%% qflow from raw (affine registered to qflow_ap)
qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_raw_affine_to_ap';
cd(qDir);

phFiles = dir('*_ph_corr_affine.nii.gz');
phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl

for ii = 1:numel(phFiles)    
    nii = load_nii(phFiles{ii});   
    PH(:,:,:,ii) = nii.img;
end

venc = 30;
VEL = venc * (PH./pi);

vq = double(VEL);

% mask
M = load_nii('../mask/inside_pipes_mask3.nii.gz');
M = single(M.img);

% implay_RR([M.*vq(:,:,:,1), ...
%            M.*vq(:,:,:,2), ...
%            M.*vq(:,:,:,3)] ...
%            ,'jet',[-5,5]);


%% qflow from raw with Rigid
% qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_raw_rigid';
% cd(qDir);
% 
% phFiles = dir('*_ph_corr_rigid.nii.gz');
% phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
% 
% for ii = 1:numel(phFiles)    
%     nii = load_nii(phFiles{ii});   
%     PH(:,:,:,ii) = nii.img;
% end
% 
% venc = 30;
% VEL = venc * (PH./pi);
% 
% vq = double(VEL);
% 
% % mask
% M = load_nii('../mask/inside_pipes_mask.nii.gz');
% M = single(M.img);
% 
% % implay_RR([M.*vq(:,:,:,1), ...
% %            M.*vq(:,:,:,2), ...
% %            M.*vq(:,:,:,3)] ...
% %            ,'jet',[-5,5]);


%% qflow from raw with FFD
% qDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\QFlow\data_raw_ffd';
% cd(qDir);
% 
% phFiles = dir('*_ph_corr_ffd.nii.gz');
% phFiles = {phFiles([1,2,3]).name}; %re-order: ap/fh/rl
% 
% for ii = 1:numel(phFiles)    
%     nii = load_nii(phFiles{ii});   
%     PH(:,:,:,ii) = nii.img;
% end
% 
% venc = 30;
% VEL = venc * (PH./pi);
% 
% vq = double(VEL);
% 
% % mask
% M = load_nii('../mask/inside_pipes_mask.nii.gz');
% M = single(M.img);
% 
% % implay_RR([M.*vq(:,:,:,1), ...
% %            M.*vq(:,:,:,2), ...
% %            M.*vq(:,:,:,3)] ...
% %            ,'jet',[-5,5]);


%% bssfp transformed only
bDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\data_trans';
cd(bDir);
% cd 3_phase_stacks_recon
% cd 5_phase_stacks_recon
cd 6_phase_stacks_recon

nii = load_nii('ap_vel_trans.nii.gz');
vb(:,:,:,1) = double(nii.img);
nii = load_nii('fh_vel_trans.nii.gz');
vb(:,:,:,2) = -1.*double(nii.img);
nii = load_nii('rl_vel_trans.nii.gz');
vb(:,:,:,3) = double(nii.img);


% masks
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\mask');
P = load_nii('inside_pipes_mask.nii.gz');
S = load_nii('spherical_part_mask.nii.gz');

M = double(P.img & S.img);

vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% bssfp rigid
% bDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\data_rigid';
% cd(bDir);
% cd 5_phase_stacks_recon
% 
% nii = load_nii('ap_vel_rigid.nii.gz');
% vb(:,:,:,1) = double(nii.img);
% nii = load_nii('fh_vel_rigid.nii.gz');
% vb(:,:,:,2) = -1.*double(nii.img);
% nii = load_nii('rl_vel_rigid.nii.gz');
% vb(:,:,:,3) = double(nii.img);
% 
% 
% % masks
% cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\mask');
% P = load_nii('inside_pipes_mask.nii.gz');
% S = load_nii('spherical_part_mask.nii.gz');
% 
% M = double(P.img & S.img);
% 
% vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
% vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% bssfp ffd
% bDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\data_ffd';
% cd(bDir);
% cd 6_phase_stacks_recon
% 
% nii = load_nii('ap_vel_ffd.nii.gz');
% vb(:,:,:,1) = double(nii.img);
% nii = load_nii('fh_vel_ffd.nii.gz');
% vb(:,:,:,2) = -1.*double(nii.img);
% nii = load_nii('rl_vel_ffd.nii.gz');
% vb(:,:,:,3) = double(nii.img);
% 
% 
% % masks
% cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\mask');
% P = load_nii('inside_pipes_mask.nii.gz');
% S = load_nii('spherical_part_mask.nii.gz');
% 
% M = double(P.img & S.img);
% 
% vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
% vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% bssfp rigid to qflow_ap
bDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\data_rigid_to_qflow';
cd(bDir);
% cd 3_phase_stacks_recon
% cd 5_phase_stacks_recon
% cd 6_phase_stacks_recon
cd 5_phase_stacks_recon_noSag

nii = load_nii('ap_vel_rigid.nii.gz');
vb(:,:,:,1) = double(nii.img);
nii = load_nii('fh_vel_rigid.nii.gz');
vb(:,:,:,2) = -1.*double(nii.img);
nii = load_nii('rl_vel_rigid.nii.gz');
vb(:,:,:,3) = double(nii.img);


% masks
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\mask');
P = load_nii('inside_pipes_mask3.nii.gz');
S = load_nii('spherical_part_mask.nii.gz');

% M = double(P.img & S.img);
M = double(P.img);

vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% bssfp affine to qflow_ap
bDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\data_affine_to_qflow';
cd(bDir);
% cd 3_phase_stacks_recon
% cd 5_phase_stacks_recon
% cd 6_phase_stacks_recon
cd 5_phase_stacks_recon_noSag
% cd 4_phase_stacks_recon_noSag

nii = load_nii('ap_vel_affine.nii.gz');
vb(:,:,:,1) = double(nii.img);
nii = load_nii('fh_vel_affine.nii.gz');
vb(:,:,:,2) = -1.*double(nii.img);
nii = load_nii('rl_vel_affine.nii.gz');
vb(:,:,:,3) = double(nii.img);


% masks
cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_21_Flow_Phantom9\PC-bSSFP\mask');
P = load_nii('inside_pipes_mask3.nii.gz');
% S = load_nii('spherical_part_mask.nii.gz');

M = double(P.img);

vq3d = sqrt(vq(:,:,:,1).^2 + vq(:,:,:,2).^2 + vq(:,:,:,3).^2);
vb3d = sqrt(vb(:,:,:,1).^2 + vb(:,:,:,2).^2 + vb(:,:,:,3).^2);


%% View velocity images
implay_RR([M.*vq(:,:,:,1), M.*vq(:,:,:,2), M.*vq(:,:,:,3), M.*vq3d;...
           M.*vb(:,:,:,1), M.*vb(:,:,:,2), M.*vb(:,:,:,3), M.*vb3d]...
           ,'jet',[-10,10]);

       
       

%% Plot histograms

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

figure;

% set background/really low values to NaN
idx1 = find(M_vq(:,:,:,1)>-bkrdClip & M_vq(:,:,:,1)<bkrdClip); %M_vq(idx) = NaN;
idx1 = find(M_vb(:,:,:,1)>-bkrdClip & M_vb(:,:,:,1)<bkrdClip); %M_vb(idx) = NaN;
idx2 = find(M_vq(:,:,:,2)>-bkrdClip & M_vq(:,:,:,2)<bkrdClip);
idx2 = find(M_vb(:,:,:,2)>-bkrdClip & M_vb(:,:,:,2)<bkrdClip);
idx3 = find(M_vq(:,:,:,3)>-bkrdClip & M_vq(:,:,:,3)<bkrdClip);
idx3 = find(M_vb(:,:,:,3)>-bkrdClip & M_vb(:,:,:,3)<bkrdClip);

idx1 = find(M_vq(:,:,:,1)>-bkrdClip & M_vq(:,:,:,1)<bkrdClip);
idx2 = find(M_vq(:,:,:,2)>-bkrdClip & M_vq(:,:,:,2)<bkrdClip);
idx3 = find(M_vq(:,:,:,3)>-bkrdClip & M_vq(:,:,:,3)<bkrdClip);
idx_unique = unique([idx1; idx2; idx3]);
M_vq(idx_unique) = NaN; nanMask_vq = isnan(M_vq(:,:,:,1));
M_vb(idx_unique) = NaN; nanMask_vb = isnan(M_vb(:,:,:,1));
% M_vq3d(idx_unique) = NaN;
% M_vb3d(idx_unique) = NaN;

% set background/really low values to NaN
% idx = find(M_vq>-bkrdClip & M_vq<bkrdClip); M_vq(idx) = NaN;
% idx = find(M_vb>-bkrdClip & M_vb<bkrdClip); M_vb(idx) = NaN;
% idx = find(M_vq3d>-bkrdClip & M_vq3d<bkrdClip); M_vq3d(idx) = NaN;
% idx = find(M_vb3d>-bkrdClip & M_vb3d<bkrdClip); M_vb3d(idx) = NaN;

% % display with thresholding
% implay_RR([M_vq(:,:,:,1), M_vq(:,:,:,2), M_vq(:,:,:,3), M_vq3d;...
%            M_vb(:,:,:,1), M_vb(:,:,:,2), M_vb(:,:,:,3), M_vb3d]...
%            ,'jet',[-10,10]);


for ii = 1:3
    subplot(2,2,ii);
    temp = M_vq(:,:,:,ii); disp(num2str(numel(temp)));
    QVhist(:,ii) = temp(~nanMask_vq); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
    histogram(QVhist(:,ii),'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');

    hold on;

    temp = M_vb(:,:,:,ii); disp(num2str(numel(temp)));
    BVhist(:,ii) = temp(~nanMask_vq); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
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

lowerThresh = 0.5;
upperThresh = 20;
binSize = 0.1;

subplot(2,2,4);
QVhist(:,4) = M_vq3d(~nanMask_vq); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
histogram(QVhist(:,4),'binEdges',0:binSize:upperThresh,'Normalization','Probability');

hold on;

BVhist(:,4) = M_vb3d(~nanMask_vq); %QBhist = QBhist(QBhist>lowerThresh); QBhist = QBhist(QBhist<upperThresh);
histogram(BVhist(:,4),'binEdges',0:binSize:upperThresh,'Normalization','Probability');

legend('PC-SPGR','PC-bSSFP');
xlabel('Velocity [cm/s]');
ylabel('Counts');
title('Magnitude');

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
BA.ap(:,1) = QVhist(~isnan(BVhist(:,1))); %ap
BA.ap(:,2) = BVhist(~isnan(BVhist(:,1)));
BA.fh(:,1) = QVhist(~isnan(BVhist(:,2))); %fh
BA.fh(:,2) = BVhist(~isnan(BVhist(:,2)));
BA.rl(:,1) = QVhist(~isnan(BVhist(:,3))); %rl
BA.rl(:,2) = BVhist(~isnan(BVhist(:,3)));
BA.mag(:,1) = QVhist(~isnan(BVhist(:,4))); %magnitude
BA.mag(:,2) = BVhist(~isnan(BVhist(:,4)));
% length(column)


%% Voxel-wise difference

lowerThresh = -10;
upperThresh = 10;
binSize = 0.1;

bkrdClip = 0.025;

for ii = 1:3
    M_vq(:,:,:,ii) = M.*vq(:,:,:,ii);
    M_vb(:,:,:,ii) = M.*vb(:,:,:,ii);
end

M_vq3d = M.*vq3d;
M_vb3d = M.*vb3d;

figure;

% set background/really low values to NaN
idx = find(M_vq>-bkrdClip & M_vq<bkrdClip); M_vq(idx) = NaN;
idx = find(M_vb>-bkrdClip & M_vb<bkrdClip); M_vb(idx) = NaN;
idx = find(M_vq3d>-bkrdClip & M_vq3d<bkrdClip); M_vq3d(idx) = NaN;
idx = find(M_vb3d>-bkrdClip & M_vb3d<bkrdClip); M_vb3d(idx) = NaN;

% % display with thresholding
% implay_RR([M_vq(:,:,:,1), M_vq(:,:,:,2), M_vq(:,:,:,3), M_vq3d;...
%            M_vb(:,:,:,1), M_vb(:,:,:,2), M_vb(:,:,:,3), M_vb3d]...
%            ,'jet',[-10,10]);


for ii = 1:3
    subplot(2,2,ii);
    temp = M_vq(:,:,:,ii)-M_vb(:,:,:,ii); disp(num2str(numel(temp)));
    temp2(:,ii)=temp(:);
    DIFFhist = temp(:); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
    histogram(DIFFhist,'binEdges',lowerThresh:binSize:upperThresh,'Normalization','Probability');

    legend('PC-SPGR - PC-bSSFP');
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

lowerThresh = 0.5;
upperThresh = 20;
binSize = 0.1;

subplot(2,2,4);
DIFFhist = M_vq3d(:)-M_vb3d(:); %QVhist = QVhist(QVhist>lowerThresh); QVhist = QVhist(QVhist<upperThresh);
temp2(:,4)=DIFFhist(:);
histogram(DIFFhist(:),'binEdges',0:binSize:upperThresh,'Normalization','Probability');

legend('PC-SPGR - PC-bSSFP');
xlabel('Velocity [cm/s]');
ylabel('Counts');
title('Magnitude');