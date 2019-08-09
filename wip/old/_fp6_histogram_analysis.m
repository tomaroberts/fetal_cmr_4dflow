%% Flow Phantom 6 - histogram analysis

%% directory containing cropped volumes
cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\#Finalised_Data\flow_phantom6_bssfp_ffd_registered2qflow\#histogram_analysis');

QV = load_nii('QFlow_VEL_ffd_Vx_cropped.nii.gz');
temp = load_nii('QFlow_VEL_ffd_Vy_cropped.nii.gz'); QV.img(:,:,:,2) = temp.img;
temp = load_nii('QFlow_VEL_ffd_Vz_cropped.nii.gz'); QV.img(:,:,:,3) = temp.img;

MPQ = load_nii('QFlow_MAG_ffd_cropped_pipe_mask.nii.gz');

BV = load_nii('vel_vol_ffd_bssfp2qflow_ap_cropped.nii.gz');
temp = load_nii('vel_vol_ffd_bssfp2qflow_fh_cropped.nii.gz'); BV.img(:,:,:,2) = temp.img;
temp = load_nii('vel_vol_ffd_bssfp2qflow_rl_cropped.nii.gz'); BV.img(:,:,:,3) = temp.img;

MF = load_nii('vel_vol_ffd_bssfp2qflow_cropped_flask_mask.nii.gz');

MPB1 = load_nii('vel_vol_ffd_bssfp2qflow_cropped_pipe_mask1.nii.gz');
MPB2 = load_nii('vel_vol_ffd_bssfp2qflow_cropped_pipe_mask2.nii.gz');
MPB3 = load_nii('vel_vol_ffd_bssfp2qflow_cropped_pipe_mask3.nii.gz');

MPB.img = double(logical(MPB1.img + MPB2.img + MPB3.img)) .* double(MF.img) ;

perm = [2,1,3,4];
QV.img = double(flip(flip(permute(QV.img,perm),1),2));
BV.img = double(flip(flip(permute(BV.img,perm),1),2));
MF.img = double(flip(flip(permute(MF.img,perm),1),2));
MPQ.img = double(flip(flip(permute(MPQ.img,perm),1),2));
MPB.img = double(flip(flip(permute(MPB.img,perm),1),2));

% rearrange bssfp components to match Qflow
BV.img = cat(4,BV.img(:,:,:,3),BV.img(:,:,:,1),BV.img(:,:,:,2));

% implay_RR([ QV.img(:,:,:,1) QV.img(:,:,:,2) QV.img(:,:,:,3); ...
%             MF.img.*BV.img(:,:,:,1) MF.img.*BV.img(:,:,:,2) MF.img.*BV.img(:,:,:,3)], 'jet', [-50,50]);

QVmag.img = sqrt( QV.img(:,:,:,1).^2 + QV.img(:,:,:,3).^2 + QV.img(:,:,:,3).^2 );
BVmag.img = sqrt( BV.img(:,:,:,1).^2 + BV.img(:,:,:,3).^2 + BV.img(:,:,:,3).^2 );

implay_RR([ MPQ.img.*QVmag.img MPB.img.*MF.img.*BVmag.img ], 'jet', [-50,50]);



%% Plot histograms
noiseThresh = 0;

figure;
QVhist = MPQ.img.*QVmag.img;
QVhist = QVhist(:); QVhist = QVhist(QVhist>noiseThresh);
histogram(QVhist,50,'Normalization','Probability');
% histogram(QVhist,50);

hold on;

% BVhist = MPQ.img.*BVmag.img;
BVhist = MPB.img.*MF.img.*BVmag.img;
BVhist = BVhist(:); BVhist = BVhist(BVhist>noiseThresh);
histogram(BVhist,50,'Normalization','Probability');
% histogram(BVhist,50);

% Qpdca = fitdist(QVhist,'Normal');
% Bpdca = fitdist(BVhist,'Normal');
% plot(0:50,pdf(Qpdca,0:50));
% plot(0:50,pdf(Bpdca,0:50));
hold on;

Qgm = fitgmdist(QVhist,2);
Qgmdist = gmdistribution(Qgm.mu,Qgm.Sigma,Qgm.ComponentProportion);
plot([0:0.25:50]',pdf(Qgmdist,[0:0.25:50]'),'b-');

Bgm = fitgmdist(BVhist,2);
Bgmdist = gmdistribution(Bgm.mu,Bgm.Sigma,Bgm.ComponentProportion);
plot([0:0.25:50]',pdf(Bgmdist,[0:0.25:50]'),'r-');

legend('PC-SPGR','PC-bSSFP', ...
       ['mu = ' num2str(min(Qgm.mu)) ' / ' num2str(max(Qgm.mu))], ...
       ['mu = ' num2str(min(Bgm.mu)) ' / ' num2str(max(Bgm.mu))]);
xlabel('Velocity [cm/s]');
ylabel('Counts');





       