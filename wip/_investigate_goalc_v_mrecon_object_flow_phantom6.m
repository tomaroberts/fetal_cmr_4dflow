%% get parameters from PC-bSSFP files
%
% Flow_Phantom_6 dataset
% - Trying to collect the parameters from the raw files that I need to
% correctly orient the slices and velocity-encodings

sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP';
cd(sp);

%% file names
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


%% load_untouch_nii:

tra.p.nii = load_untouch_nii([scanName{1} '_MAG.nii.gz']);
sag.p.nii = load_untouch_nii([scanName{2} '_MAG.nii.gz']);
cor.p.nii = load_untouch_nii([scanName{3} '_MAG.nii.gz']);
tra.m.nii = load_untouch_nii([scanName{4} '_MAG.nii.gz']);
sag.m.nii = load_untouch_nii([scanName{5} '_MAG.nii.gz']);
cor.m.nii = load_untouch_nii([scanName{6} '_MAG.nii.gz']);
tra.o.nii = load_untouch_nii([scanName{7} '_MAG.nii.gz']);
sag.o.nii = load_untouch_nii([scanName{8} '_MAG.nii.gz']);
cor.o.nii = load_untouch_nii([scanName{9} '_MAG.nii.gz']);



%% get values from goalc.txt files

% GR
varArray = { 'ori_x', 'ori_y', 'ori_z', ...          % variables to collect from goalc files
             'flow1', 'flow2'};

objects = { 'GR `m[0]', 'GR `mc[0]', 'GR `m[3]', ... % SLICE-SELECT, S, gradients
            'GR `s_ex', 'GR `r[0]', ...              % READOUT, M, gradients
            };                 

objStructNames = { 'm0', 'mc', 'm3', ...             % cell array for creating structures
                   's_ex', 'r0' };

% get variables for each scan
for ii = 1:numel(objects)

    tra.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{1} '_goalc.txt'], objects{ii}, varArray );
    sag.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{2} '_goalc.txt'], objects{ii}, varArray );
    cor.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{3} '_goalc.txt'], objects{ii}, varArray );

    tra.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{4} '_goalc.txt'], objects{ii}, varArray );
    sag.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{5} '_goalc.txt'], objects{ii}, varArray );
    cor.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{6} '_goalc.txt'], objects{ii}, varArray );

    tra.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{7} '_goalc.txt'], objects{ii}, varArray );
    sag.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{8} '_goalc.txt'], objects{ii}, varArray );
    cor.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{9} '_goalc.txt'], objects{ii}, varArray );

end


% LCA
varArray = { 'orientation:[0]', 'prep_dir:[0]', 'fat_shift_dir:[0]', 'pat_angulations:[0]' }; % variables to collect from goalc files

objects = { 'LCA `ima' }; % LCA object                 

objStructNames = { 'LCA' }; % cell array for creating structures

% get variables for each scan
for ii = 1:numel(objects)

    tra.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{1} '_goalc.txt'], objects{ii}, varArray );
    sag.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{2} '_goalc.txt'], objects{ii}, varArray );
    cor.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{3} '_goalc.txt'], objects{ii}, varArray );

    tra.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{4} '_goalc.txt'], objects{ii}, varArray );
    sag.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{5} '_goalc.txt'], objects{ii}, varArray );
    cor.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{6} '_goalc.txt'], objects{ii}, varArray );

    tra.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{7} '_goalc.txt'], objects{ii}, varArray );
    sag.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{8} '_goalc.txt'], objects{ii}, varArray );
    cor.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{9} '_goalc.txt'], objects{ii}, varArray );

end


% O_MATRIX
varArray = { 'm_orient', 'p_orient', 's_orient' }; % variables to collect from goalc files

objects = { 'O_MATRIX `locations[0]' }; % LCA object                 

objStructNames = { 'O_MATRIX_loc0' }; % cell array for creating structures

% get variables for each scan
for ii = 1:numel(objects)

    tra.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{1} '_goalc.txt'], objects{ii}, varArray );
    sag.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{2} '_goalc.txt'], objects{ii}, varArray );
    cor.p.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{3} '_goalc.txt'], objects{ii}, varArray );

    tra.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{4} '_goalc.txt'], objects{ii}, varArray );
    sag.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{5} '_goalc.txt'], objects{ii}, varArray );
    cor.m.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{6} '_goalc.txt'], objects{ii}, varArray );

    tra.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{7} '_goalc.txt'], objects{ii}, varArray );
    sag.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{8} '_goalc.txt'], objects{ii}, varArray );
    cor.o.HDR.(objStructNames{ii}) = get_goalc_parameters( [scanName{9} '_goalc.txt'], objects{ii}, varArray );

end


%% get relevant MRecon parameters

[tra.p.PARAM, tra.p.A, tra.p.Anorm] = get_mrecon_params(scanName{1});
[sag.p.PARAM, sag.p.A, sag.p.Anorm] = get_mrecon_params(scanName{2});
[cor.p.PARAM, cor.p.A, cor.p.Anorm] = get_mrecon_params(scanName{3});

[tra.m.PARAM, tra.m.A, tra.m.Anorm] = get_mrecon_params(scanName{4});
[sag.m.PARAM, sag.m.A, sag.m.Anorm] = get_mrecon_params(scanName{5});
[cor.m.PARAM, cor.m.A, cor.m.Anorm] = get_mrecon_params(scanName{6});

[tra.o.PARAM, tra.o.A, tra.o.Anorm] = get_mrecon_params(scanName{7});
[sag.o.PARAM, sag.o.A, sag.o.Anorm] = get_mrecon_params(scanName{8});
[cor.o.PARAM, cor.o.A, cor.o.Anorm] = get_mrecon_params(scanName{9});




%% subfn to get params from PARAM.mat and the .nii affine
function [P,A,Anorm] = get_mrecon_params(scanName)

    load([scanName '_PARAM.mat']);

    P.Angulation = PARAM.Scan.Angulation;
    P.Orientation = PARAM.Scan.Orientation;
    P.FoldOverDir = PARAM.Scan.FoldOverDir;
    P.FatShiftDir = PARAM.Scan.FatShiftDir;
    P.ijk = PARAM.Scan.ijk;
    P.MPS = PARAM.Scan.MPS;
    P.xyz = PARAM.Scan.xyz;
    P.REC = PARAM.Scan.REC;

    nii = load_untouch_nii([scanName '_MAG.nii.gz']);

    P.pixdim = nii.hdr.dime.pixdim(2:4);
    P.dim    = nii.hdr.dime.dim(2:4);
    P.FOV    = P.pixdim .* P.dim;
    
    A(1,:) = nii.hdr.hist.srow_x;
    A(2,:) = nii.hdr.hist.srow_y;
    A(3,:) = nii.hdr.hist.srow_z;
    A(4,:) = [0 0 0 1];
    
    Anorm = A;
    Anorm(:,1) = Anorm(:,1) ./ P.pixdim(1);
    Anorm(:,2) = Anorm(:,2) ./ P.pixdim(2);
    Anorm(:,3) = Anorm(:,3) ./ P.pixdim(3);
    Anorm(:,4) = [];
    Anorm(4,:) = [];

end






