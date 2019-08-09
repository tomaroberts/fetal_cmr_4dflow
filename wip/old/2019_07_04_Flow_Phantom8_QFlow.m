%% 2019_07_04 --- QFlow --- FlowPhantom_6 --- spherical flask
%
% nb: recon from raw because Philips phase processing too severe in .nii
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1)
%% Files reconstructed on beastie02 using FP8_pc_mrecon_from_raw.m
%- outputs: _RE / _IM / _CPLX / _MAG / _PH.nii.gz files + others
reconDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_04_Flow_Phantom8\QFlow\data';
cd(reconDir);

%% Step 2)
%% Make directory stucture and copy magnitude/complex data
% In bash:
% mkdir data
% mkdir mask
%
% cp *_MAG* ../data/
% cp *_CPLX* ../data/
% cp *_RE* ../data/
% cp *_IM* ../data/
%
% Rename accordingly...
%
% Make ../mask/static_water_mask.nii.gz in MITK
%

%% Step 3)
%% Make phase data
cd(reconDir);
fileNames = dir('*_re_trans_qflow2bssfp128px.nii.gz');

for ii = 1:numel(fileNames)
    currStr = erase(fileNames(ii).name, 're_trans_qflow2bssfp128px.nii.gz');
    
    re = load_untouch_nii( [currStr 're_trans_qflow2bssfp128px.nii.gz'] );
    im = load_untouch_nii( [currStr 'im_trans_qflow2bssfp128px.nii.gz'] );
    
    % initialise ph/ph_corr nifti structures
    ph = re; ph_corr = re;
        
    % Create uncorrected phase images
    cx.img = re.img + 1i.*im.img;
    ph.img = angle(abs(cx.img).*exp(sqrt(-1)*(angle(cx.img)+pi))); %DO I NEED +pi in exp? Not sure

    % save ph (without polynomial correction)
    ph.img = single(ph.img); % convert to single to match original .nii
    save_untouch_nii(ph,[currStr 'ph_trans_qflow2bssfp128px.nii.gz']);

    % Create polynomial corrected phase images
    cd ../../PC-bSSFP/mask/
    static_water_mask = load_untouch_nii('static_water_mask_128px.nii.gz');

    % run polynomial correction
    % - Y is corrected complex data, ie: cx_corr.img
    [Y, P0, P1] = phase_correction_poly_phantom( cx.img, ...
                               'staticwatermask', logical(static_water_mask.img) );
    ph_corr.img = angle(abs(Y).*exp(sqrt(-1)*(angle(Y))));

    % save ph_corr.nii.gz
    ph_corr.img = single(ph_corr.img); % convert to single to match original .nii
    
    cd ../../QFlow/data/
    save_untouch_nii(ph_corr,[currStr 'ph_corr_trans_qflow2bssfp128px.nii.gz']);
end

%% Step 4)
%% Load stacks
reconDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\2019_07_04_Flow_Phantom8\QFlow';
cd(reconDir);
cd data

magFiles = dir('*_mag_*.nii.gz');
magFiles = {magFiles([1,2,3]).name}; %re-order: ap/fh/rl

phFiles = dir('*_ph_corr_*.nii.gz');
phFiles = {phFiles([1,2,3]).name};

for ii = 1:numel(magFiles)
    nii = load_untouch_nii(magFiles{ii});
    MAG(:,:,:,ii) = nii.img;
    
    nii = load_untouch_nii(phFiles{ii});
    PH(:,:,:,ii) = nii.img;
end

clear nii

% mask
M = load_untouch_nii('../mask/inside_pipes_mask_trans_qflow2bssfp128px.nii.gz');
M = single(M.img);

% implay_RR([M.*PH(:,:,:,1), M.*PH(:,:,:,2), M.*PH(:,:,:,3)]...
%            ,'gray',[-pi/6,pi/6]);

venc = 30;
VEL = venc * (PH./pi);
       
%% Step 4)
%% View velocity volumes  
permOrder = [2,1,3];

% implay_RR([permute(MAG(:,:,:,1),permOrder), ...
%            permute(MAG(:,:,:,2),permOrder), ...
%            permute(MAG(:,:,:,3),permOrder)] ...
%            ,'gray');

implay_RR([permute(M.*PH(:,:,:,1),permOrder), ...
           permute(M.*PH(:,:,:,2),permOrder), ...
           permute(M.*PH(:,:,:,3),permOrder)] ...
           ,'jet',[-pi/6,pi/6]);

implay_RR([permute(M.*VEL(:,:,:,1),permOrder), ...
           permute(M.*VEL(:,:,:,2),permOrder), ...
           permute(M.*VEL(:,:,:,3),permOrder)] ...
           ,'jet',[-5,5]);
       
