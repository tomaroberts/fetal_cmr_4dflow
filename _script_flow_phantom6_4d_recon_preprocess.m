%% rename files for volume reconstruction

sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP\raw_data';
cd(sp);

mkdir ../data
mkdir ../mask
mkdir ../flow
mkdir ../ph_vol_rl
mkdir ../ph_vol_ap
mkdir ../ph_vol_fh


%% filenames

% M+
scanName{1} = 'fl_01112018_2023463_30_2_pih1bffem2dmplusaxclearV4';
scanName{2} = 'fl_01112018_2029121_31_2_pih1bffem2dmplussagclearV4';
scanName{3} = 'fl_01112018_2031445_32_2_pih1bffem2dmpluscorclearV4';

% M-
scanName{4} = 'fl_01112018_2057445_36_2_pih1bffem2dmminusaxexperiment1clearV4';
scanName{5} = 'fl_01112018_2102065_37_2_pih1bffem2dmminussagexperiment1clearV4';
scanName{6} = 'fl_01112018_2105508_38_2_pih1bffem2dmminuscorexperiment1clearV4';

% Oblique
scanName{7} = 'fl_01112018_2042459_33_2_pih1bffem2dmplusaxrotap45clearV4';
scanName{8} = 'fl_01112018_2048374_34_2_pih1bffem2dmplussagrotrl45clearV4';
scanName{9} = 'fl_01112018_2050360_35_2_pih1bffem2dmpluscorrotap45clearV4';

newNames = {'tp','sp','cp',...
            'tm','sm','cm',...
            'to','so','co'};

% Static Water Mask
MASK = load_untouch_nii('../mask/static_water_mask.nii.gz');
MASK = single(MASK.img);
              
              
%% rename, move to data folder

for ii = 1:numel(scanName)
    copyfile([scanName{ii} '_MAG.nii.gz'], ['../data/' newNames{ii} '_mag.nii.gz']);
    copyfile([scanName{ii} '_CPLX.nii.gz'],['../data/' newNames{ii} '_cplx.nii.gz']);
end


%% make phase data, save to data folder

for ii = 1:numel(scanName)
    
    % apply polynomial phase correction
    nii = load_untouch_nii([scanName{ii} '_PH.nii.gz']);
    nii.img = poly3D_phase_subtract(nii.img,MASK);
    save_untouch_nii(nii,['../data/' newNames{ii} '_ph.nii.gz']);
    
    % offset images by pi for 4d recon pipeline
    nii.img = nii.img + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset.nii.gz']);
end


%% make component phase images in world coordinates

% First moment values from GVE in MPS coordinates
Vm = 9.95;
Vp = 0;
Vs = -10.98;

Vmps = [Vm; Vp; Vs];


for ii = 1:numel(scanName)
        
    % load existing PH image, replace img, keep header same (to maintain affine matrix)
    nii = load_untouch_nii([scanName{ii} '_PH.nii.gz']);
    nii.img = poly3D_phase_subtract(nii.img,MASK);
    
    % get Vxyz/Vworld (currently don't use these)
    [Vworld(:,ii), Vxyz(:,ii)] = vmps2vworld(Vmps,get_pcmr_orientation_parameters( [scanName{ii} '_goalc.txt'] ));
        
    % separate PH image into components
    [Pcomp, Vworld_frac(ii,:), Vxyz_frac(ii,:)] = phase2components(nii.img,Vmps,[scanName{ii} '_goalc.txt']);
    
    % Pcomp - dim5 = 1 => Vworld, dim5 = 2 => Vxyz
    % save Vworld_component images (+ offset)
    nii.img = Pcomp(:,:,:,1,1) + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_rl.nii.gz']);
    nii.img = Pcomp(:,:,:,2,1) + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_ap.nii.gz']);
    nii.img = Pcomp(:,:,:,3,1) + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_fh.nii.gz']);
    
    % save Vxyz_component images (+ offset)
    nii.img = Pcomp(:,:,:,1,2) + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_x.nii.gz']);
    nii.img = Pcomp(:,:,:,2,2) + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_y.nii.gz']);
    nii.img = Pcomp(:,:,:,3,2) + pi;
    save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_z.nii.gz']);

end

disp('Finished making component phase images.');


%% create frameweight .txt files
%-tra.p/sag.p/cor.p/tra.o/sag.o
weights = abs(Vworld_frac([1,2,3,7,8],:));
weights(:,:,2) = abs(Vxyz_frac([1,2,3,7,8],:));
nS = size(nii.img,3);

% RL
weights_rl = repmat(weights(:,1,1),1,nS);
weights_rl = reshape(weights_rl',1,numel(weights_rl));
dlmwrite( '../flow/frame_weight_rl.txt', weights_rl, 'delimiter', ' ' );
weights_x = repmat(weights(:,1,2),1,nS);
weights_x = reshape(weights_x',1,numel(weights_x));
dlmwrite( '../flow/frame_weight_x.txt', weights_x, 'delimiter', ' ' );

% AP
weights_ap = repmat(weights(:,2,1),1,nS);
weights_ap = reshape(weights_ap',1,numel(weights_ap));
dlmwrite( '../flow/frame_weight_ap.txt', weights_ap, 'delimiter', ' ' );
weights_y = repmat(weights(:,2,2),1,nS);
weights_y = reshape(weights_y',1,numel(weights_y));
dlmwrite( '../flow/frame_weight_y.txt', weights_y, 'delimiter', ' ' );

% FH
weights_fh = repmat(weights(:,3,1),1,nS);
weights_fh = reshape(weights_fh',1,numel(weights_fh));
dlmwrite( '../flow/frame_weight_fh.txt', weights_fh, 'delimiter', ' ' );
weights_z = repmat(weights(:,3,2),1,nS);
weights_z = reshape(weights_z',1,numel(weights_z));
dlmwrite( '../flow/frame_weight_z.txt', weights_z, 'delimiter', ' ' );


%% SEND THROUGH VOLUME RECON:
%- get out: 3 phase component volumes
cd(sp);
copyfile('../ph_vol_rl/ph_vol_rl.nii.gz','../flow');
copyfile('../ph_vol_ap/ph_vol_ap.nii.gz','../flow');
copyfile('../ph_vol_fh/ph_vol_fh.nii.gz','../flow');


%% convert rl/ap/fh-component volumes to 3D.nii
cd(sp);
cd ../data

nii_mag_vol   = load_untouch_nii('../mag_vol/mag_vol.nii.gz');
nii_ph_vol_rl = load_untouch_nii('../flow/ph_vol_rl.nii.gz');
nii_ph_vol_ap = load_untouch_nii('../flow/ph_vol_ap.nii.gz');
nii_ph_vol_fh = load_untouch_nii('../flow/ph_vol_fh.nii.gz');

% create new .nii based on one of the existing ph_vol
%%% RL-AP-FH
nii_ph_vol = nii_ph_vol_rl;
nii_ph_vol.img(:,:,:,2) = nii_ph_vol_ap.img;
nii_ph_vol.img(:,:,:,3) = nii_ph_vol_fh.img;
extWorld = 'rl_ap_hf';

%%% RL-FH-AP
% nii_ph_vol = nii_ph_vol_rl;
% nii_ph_vol.img(:,:,:,2) = nii_ph_vol_fh.img;
% nii_ph_vol.img(:,:,:,3) = nii_ph_vol_ap.img;
% extWorld = 'rl_fh_ap';

%%% AP-RL-FH
% nii_ph_vol = nii_ph_vol_ap;
% nii_ph_vol.img(:,:,:,2) = nii_ph_vol_rl.img;
% nii_ph_vol.img(:,:,:,3) = nii_ph_vol_fh.img;
% extWorld = 'ap_rl_fh';

%%% AP-FH-RL
% nii_ph_vol = nii_ph_vol_ap;
% nii_ph_vol.img(:,:,:,2) = nii_ph_vol_fh.img;
% nii_ph_vol.img(:,:,:,3) = nii_ph_vol_rl.img;
% extWorld = 'ap_fh_rl';

%%% FH-RL-AP
% nii_ph_vol = nii_ph_vol_fh;
% nii_ph_vol.img(:,:,:,2) = nii_ph_vol_rl.img;
% nii_ph_vol.img(:,:,:,3) = nii_ph_vol_ap.img;
% extWorld = 'fh_rl_ap';

%%% FH-AP-RL
% nii_ph_vol = nii_ph_vol_fh;
% nii_ph_vol.img(:,:,:,2) = nii_ph_vol_ap.img;
% nii_ph_vol.img(:,:,:,3) = nii_ph_vol_rl.img;
% extWorld = 'fh_ap_rl';

% reset volume reconstruction phase offset
nii_ph_vol.img = nii_ph_vol.img - pi;

% adjust hdr
nii_ph_vol.hdr.dime.dim(5) = 3;
% nii_ph_vol.hdr.dime.glmax = pi; %<- unnecessary as calculated by save_untouch_nii
% nii_ph_vol.hdr.dime.glmin = -pi;

% save
save_untouch_nii(nii_ph_vol,['../flow/ph_vol_3d_' extWorld '.nii.gz']);



%% convert x/y/z-component volumes to 3D.nii
cd(sp);
cd ../data

nii_mag_vol = load_untouch_nii('../mag_vol/mag_vol.nii.gz');
nii_ph_vol_x = load_untouch_nii('../flow/ph_vol_x.nii.gz');
nii_ph_vol_y = load_untouch_nii('../flow/ph_vol_y.nii.gz');
nii_ph_vol_z = load_untouch_nii('../flow/ph_vol_z.nii.gz');

% create new .nii based on one of the existing ph_vol
%%% x-y-z
nii_ph_vol = nii_ph_vol_x;
nii_ph_vol.img(:,:,:,2) = nii_ph_vol_y.img;
nii_ph_vol.img(:,:,:,3) = nii_ph_vol_z.img;
extWorld = 'x_y_z';

% reset volume reconstruction phase offset
nii_ph_vol.img = nii_ph_vol.img - pi;

% adjust hdr
nii_ph_vol.hdr.dime.dim(5) = 3;

% save
save_untouch_nii(nii_ph_vol,['../flow/ph_vol_3d_' extWorld '.nii.gz']);