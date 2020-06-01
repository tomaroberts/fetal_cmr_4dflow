%% fcmr_4d_kt_sweep_recon
%
% Working through 4d recon pipeline with sweep data:
%
% For now: Adult dataset
%

reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\Adult\2019_03_28_Adult_Heart';
cd(reconDir);


%% Step 0)

% Standard k-t recon using MRecon
% TODO: is this working properly? Surprised it seems to work without modification


%% Convert sweep ktrecon output to '3.5D' (i.e.: time/space in 3rd dimension) volume

% TODO: fold into preproc.m ?

ktreconDir  = fullfile( reconDir, 'ktrecon' );
rltFileList = dir( fullfile( ktreconDir, 's*_rlt_ab.nii.gz' ) );
nStack = numel(rltFileList);

for iStk = 1:nStack

    R = load_untouch_nii( fullfile( ktreconDir, rltFileList(iStk).name ) );
    
    nX = size( R.img, 1 );
    nY = size( R.img, 2 );
    n3 = size( R.img, 3 ) * size( R.img, 4 );
    
    SWP.img = reshape( permute(R.img,[1,2,4,3]), nX, nY, n3 );
    
    %TODO: save new swp *.nii equivalent to LJ python output
    %TODO: retain *_rlt_* naming convention or rename *_swp_* ?
    
end

% TODO:
% copyfile( '/ktrecon', '/data' );


%% 4) Run preproc.m
cd(reconDir);
acqMethod = 'swp';
S = preproc( reconDir, acqMethod );
save( fullfile( reconDir, 'data', 'results.mat' ), 'S', '-v7.3' );
disp('Preproc complete ...');


















