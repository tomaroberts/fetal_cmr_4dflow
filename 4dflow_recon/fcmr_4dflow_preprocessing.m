function fcmr_4dflow_preprocessing( reconDir, rawDir )
%FCMR_4DFLOW_PREPROCESSING  pre-processing of 4D flow real-time data
%
%   FCMR_4DFLOW_PREPROCESSING( reconDir, ... )
%       Pre-processing required before reconstructing data with SVRTK
%       - Creates polynomial phase corrected stacks based on uterus mask
%       - Creates new goalc.txt files containing gradient first moment
%         information
% 
%   Requires:
%       Uterus ROIs drawn for each s*_.nii.gz stack, saved in ../mask/ as
%       s*_mask_uterus.nii.gz
%
%   Input:
%       reconDir            - str - fetal reconstruction directory
%
%   Optional Inputs:
%       rawDir              - str - directory on server containing .raw data
%
%   Example usage:
%       fcmr_4dflow_preprocessing( 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr194', 'Z:\' )
%
%
%   See also:
%       phase_correction_poly.m

% Tom Roberts (t.roberts@kcl.ac.uk)


%% Default paths
% - for example:

% ingenia-raw
if nargin < 2
    rawDir = 'Z:\';
end


%% Manually draw uterus ROIs for each stack for polynomial phase correction
%- save in ../mask/ as s*_mask_uterus.nii.gz

% TODO: add in check to see if mask files exist


%% Create phase (_ph.nii.gz) and phase corrected (_ph_corr.nii.gz) stacks

cd(reconDir);    
cd data
sIDs = dir('s*_rlt_ab.nii.gz');
numStacks = numel(sIDs);

cd ..    
cd ktrecon
    
for ss = 1:numStacks


    % load real/imaginary data
    re = load_untouch_nii( [sIDs(ss).name(1:3) '_rlt_re.nii.gz'] );
    im = load_untouch_nii( [sIDs(ss).name(1:3) '_rlt_im.nii.gz'] );

    % initialise ph/ph_corr nifti structures based on existing real data
    ph = re;
    ph_corr = re;


    %% Create uncorrected phase images
    cx.img = re.img + 1i.*im.img;
    ph.img = angle(abs(cx.img).*exp(sqrt(-1)*(angle(cx.img)+pi))); %DO I NEED +pi in exp? Not sure

    % save ph nifti (without polynomial correction)
    ph.img = single(ph.img); % convert to single to match original .nii
    save_untouch_nii(ph,[sIDs(ss).name(1:3) '_rlt_ph.nii.gz']);

    movefile([sIDs(ss).name(1:3) '_rlt_ph.nii.gz'] , ['../data/' sIDs(ss).name(1:3) '_rlt_ph.nii.gz']);

    disp(['Created ' sIDs(ss).name(1:3) '_rlt_ph.nii.gz']);


    %% Create polynomial corrected phase images

    cd ../mask
    uterus_mask = load_untouch_nii([sIDs(ss).name(1:3) '_mask_uterus.nii.gz']);
    heart_mask  = load_untouch_nii([sIDs(ss).name(1:3) '_mask_heart.nii.gz']);
    cd ../ktrecon

    % run polynomial correction
    % - Y is corrected complex data, ie: cx_corr.img
    [Y, P0, P1] = phase_correction_poly( cx.img, ...
                               'uterusmask', logical(uterus_mask.img), ...
                               'heartmask',  logical(heart_mask.img) );

    ph_corr.img = angle(abs(Y).*exp(sqrt(-1)*(angle(Y))));

    % save ph_corr.nii.gz
    ph_corr.img = single(ph_corr.img); % convert to single to match original .nii
    save_untouch_nii(ph_corr,[sIDs(ss).name(1:3) '_rlt_ph_corr_uterus.nii.gz']);

    movefile([sIDs(ss).name(1:3) '_rlt_ph_corr_uterus.nii.gz'] , ['../data/' sIDs(ss).name(1:3) '_rlt_ph_corr_uterus.nii.gz']);

    % save polynomial and offset
    % TODO: polynomials currently saving as 5D .nii - should be 3D
    poly_nii = re;
    polyNom = exp( -( 1i*(P0+P1) ) );

    poly_nii.img = angle(abs(polyNom).*exp(sqrt(-1)*(angle(polyNom)+pi)));

    save_untouch_nii(poly_nii,[sIDs(ss).name(1:3) '_rlt_ph_polynomial_uterus.nii.gz']);        
    movefile([sIDs(ss).name(1:3) '_rlt_ph_polynomial_uterus.nii.gz'] , ['../data/' sIDs(ss).name(1:3) '_rlt_ph_polynomial_uterus.nii.gz']); 

    disp(['Created ' sIDs(ss).name(1:3) '_rlt_ph_corr_uterus.nii.gz']);


end
    

%% Update GOALC.txt files with ORIENT information
% TODO: update ktrecon scripts so that ORIENT information included in goalc.txt file. 


end % fcmr_4dflow_preprocessing(...)