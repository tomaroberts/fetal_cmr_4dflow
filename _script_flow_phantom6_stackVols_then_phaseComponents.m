%% studypath
sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP\stack_ph_vol';
cd(sp);


%% Admin
%%% Specify which volumes to use for final reconstructions:
vols4recon = [1,2,3,7,8];

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


%% load volume reconstructed stacks
%-nb: these have already had:
%- background subtraction
%- offset applied - all in range 0-2pi.

%plus
nii.tp = load_untouch_nii(['ph_vol_' newNames{1} '.nii.gz']);
nii.sp = load_untouch_nii(['ph_vol_' newNames{2} '.nii.gz']);
nii.cp = load_untouch_nii(['ph_vol_' newNames{3} '.nii.gz']);

%minus
nii.tm = load_untouch_nii(['ph_vol_' newNames{4} '.nii.gz']);
nii.sm = load_untouch_nii(['ph_vol_' newNames{5} '.nii.gz']);
nii.cm = load_untouch_nii(['ph_vol_' newNames{6} '.nii.gz']);

%obl
nii.to = load_untouch_nii(['ph_vol_' newNames{7} '.nii.gz']);
nii.so = load_untouch_nii(['ph_vol_' newNames{8} '.nii.gz']);
nii.co = load_untouch_nii(['ph_vol_' newNames{9} '.nii.gz']);





%% Work out components of gradient moments

% First moment values from GVE in MPS coordinates
Vm = 9.95;
Vp = 0;
Vs = -10.98;

Vmps = [Vm; Vp; Vs];

Pcomp = zeros(size(nii.tp.img,1),size(nii.tp.img,2),size(nii.tp.img,3),3,numel(scanName));

for ii = 1:numel(scanName)
        
    % load existing PH vol to separate into components
    niiCurrent = load_untouch_nii(['ph_vol_' newNames{ii} '.nii.gz']);
    niiCurrent.img = niiCurrent.img - pi;
    
    % get Vxyz/Vworld (currently don't use these)
    [Vworld(:,ii), Vxyz(:,ii)] = vmps2vworld(Vmps,get_pcmr_orientation_parameters( ['../raw_data/' scanName{ii} '_goalc.txt'] ));
        
    % separate PH image into components
    [Ptemp, Vworld_frac(ii,:)] = phase2components(niiCurrent.img,Vmps,['../raw_data/' scanName{ii} '_goalc.txt']);
    
    % Pcomp - dim5 = 1 => Vworld, dim5 = 2 => Vxyz
    niiCurrent.img = Ptemp(:,:,:,1,1) + pi;
    save_untouch_nii(niiCurrent,['ph_vol_' newNames{ii} '_rl.nii.gz']);
    niiCurrent.img = Ptemp(:,:,:,2,1) + pi;
    save_untouch_nii(niiCurrent,['ph_vol_' newNames{ii} '_ap.nii.gz']);
    niiCurrent.img = Ptemp(:,:,:,3,1) + pi;
    save_untouch_nii(niiCurrent,['ph_vol_' newNames{ii} '_fh.nii.gz']);
    
    Pcomp(:,:,:,1,ii) = Ptemp(:,:,:,1,1); 
    Pcomp(:,:,:,2,ii) = Ptemp(:,:,:,2,1); 
    Pcomp(:,:,:,3,ii) = Ptemp(:,:,:,3,1);
    
%     % save Vxyz_component images (+ offset)
%     nii.img = Pcomp(:,:,:,1,2) + pi;
%     save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_x.nii.gz']);
%     nii.img = Pcomp(:,:,:,2,2) + pi;
%     save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_y.nii.gz']);
%     nii.img = Pcomp(:,:,:,3,2) + pi;
%     save_untouch_nii(nii,['../flow/' newNames{ii} '_ph_offset_z.nii.gz']);

    clear niiCurrent

end

disp('Finished making component phase images.');


%% View component phase images

cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP');

% affine:
load('Pcomp_affine_tp_sp_cp_to_so.mat');

PcompA = Pcomp;

xcrop = 109:180;
ycrop = 134:197;
s = 45;
cscale = [-pi/4,pi/4];

% cols: RL-AP-FH
imtar([Pcomp(xcrop,ycrop,s,1,1), Pcomp(xcrop,ycrop,s,2,1), Pcomp(xcrop,ycrop,s,3,1);...
       Pcomp(xcrop,ycrop,s,1,2), Pcomp(xcrop,ycrop,s,2,2), Pcomp(xcrop,ycrop,s,3,2);...
       Pcomp(xcrop,ycrop,s,1,3), Pcomp(xcrop,ycrop,s,2,3), Pcomp(xcrop,ycrop,s,3,3);...
       Pcomp(xcrop,ycrop,s,1,4), Pcomp(xcrop,ycrop,s,2,4), Pcomp(xcrop,ycrop,s,3,4);...
       Pcomp(xcrop,ycrop,s,1,5), Pcomp(xcrop,ycrop,s,2,5), Pcomp(xcrop,ycrop,s,3,5);...
       Pcomp(xcrop,ycrop,s,1,6), Pcomp(xcrop,ycrop,s,2,6), Pcomp(xcrop,ycrop,s,3,6);...
       Pcomp(xcrop,ycrop,s,1,7), Pcomp(xcrop,ycrop,s,2,7), Pcomp(xcrop,ycrop,s,3,7);...
       Pcomp(xcrop,ycrop,s,1,8), Pcomp(xcrop,ycrop,s,2,8), Pcomp(xcrop,ycrop,s,3,8);...
       Pcomp(xcrop,ycrop,s,1,9), Pcomp(xcrop,ycrop,s,2,9), Pcomp(xcrop,ycrop,s,3,9);...
       ],...
       cscale(1),cscale(2));
   

% vol_rec:
load('Pcomp_stackVols_tp_sp_cp_to_so.mat');
Pcomp = permute(Pcomp,[3,2,1,4,5]);

PcompV = Pcomp;

xcrop = 1:size(Pcomp,1);
ycrop = 1:size(Pcomp,2);
xcrop = 20:120;
ycrop = 50:136;
s = 96;
cscale = [-pi/4,pi/4];

% cols: RL-AP-FH
imtar([Pcomp(xcrop,ycrop,s,1,1), Pcomp(xcrop,ycrop,s,2,1), Pcomp(xcrop,ycrop,s,3,1);...
       Pcomp(xcrop,ycrop,s,1,2), Pcomp(xcrop,ycrop,s,2,2), Pcomp(xcrop,ycrop,s,3,2);...
       Pcomp(xcrop,ycrop,s,1,3), Pcomp(xcrop,ycrop,s,2,3), Pcomp(xcrop,ycrop,s,3,3);...
       Pcomp(xcrop,ycrop,s,1,4), Pcomp(xcrop,ycrop,s,2,4), Pcomp(xcrop,ycrop,s,3,4);...
       Pcomp(xcrop,ycrop,s,1,5), Pcomp(xcrop,ycrop,s,2,5), Pcomp(xcrop,ycrop,s,3,5);...
       Pcomp(xcrop,ycrop,s,1,6), Pcomp(xcrop,ycrop,s,2,6), Pcomp(xcrop,ycrop,s,3,6);...
       Pcomp(xcrop,ycrop,s,1,7), Pcomp(xcrop,ycrop,s,2,7), Pcomp(xcrop,ycrop,s,3,7);...
       Pcomp(xcrop,ycrop,s,1,8), Pcomp(xcrop,ycrop,s,2,8), Pcomp(xcrop,ycrop,s,3,8);...
       Pcomp(xcrop,ycrop,s,1,9), Pcomp(xcrop,ycrop,s,2,9), Pcomp(xcrop,ycrop,s,3,9);...
       ],...
       cscale(1),cscale(2));



   
% save('..\Pcomp_stackVols_tp_sp_cp_to_so','Pcomp');
   
% tom = permute(Pcomp(:,:,:,3,1),[3,2,1,4,5]);
% tom = Pcomp(:,:,:,3,1);
% implay_RR(tom);
