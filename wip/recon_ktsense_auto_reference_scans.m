%% Calculate Coil Sensitivity Maps from Dynamic Undersampled kt data

%% Admin
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr254_tom_recon';
ktDir = '\ktrecon';
sliceNum = 20;
ktFactor = 8;


%% Anon Fns
dimX = 1;
dimY = 2;

kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, dimX ), dimY ) );

phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) );
 
inv_phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [-1;+1], size(k,1)/2, 1 ) ) ) );


%% Load data
cd([reconDir ktDir]);
load(['s' num2str(sliceNum) '_csm.mat']);
load(['s' num2str(sliceNum) '_kspace.mat'],'ktAcq');


%% View mRecon-derived Coil Sensitivity Maps
sliceNum = 1;
montage_RR(squeeze(csm(:,:,1,:,sliceNum)),['Coil Sens. Maps, slice' num2str(sliceNum)],[0,1]); colormap('jet');
montage_RR(abs(imCoil(:,:,:,:,sliceNum)),['Array Coils, slice' num2str(sliceNum)],[0 5e4]);
montage_RR(abs(imBody(:,:,:,:,1)),['Body Coil, slice' num2str(sliceNum)],[0 5e4]);


%% Create Single Shot k-spaces
dim.x    = size(ktAcq,1);
dim.y    = size(ktAcq,2);
dim.dyn  = size(ktAcq,3);
dim.chan = size(ktAcq,4);
dim.loca = size(ktAcq,5);

ktAcqSS = reshape(ktAcq , dim.x, dim.y, ktFactor, dim.dyn/ktFactor, dim.chan, dim.loca );

% % view sequential k-space and associated image
% ktTom = sum(ktAcqSS(:,:,:,1,1,sliceNum),3);
% ktTom = phase_correct( ktTom );
% xtTom = kt2xt( ktTom );
% imtar(abs(ktTom),0,50);
% imtar(abs(xtTom));

ktAcqSS = phase_correct( ktAcqSS );
xtAcqSS = kt2xt( ktAcqSS );

shotNum = 1; % 1:nDyn/ktFactor

% View Single Shot Data from All Coils
montage_RR(abs(squeeze(sum(xtAcqSS(:,:,:,shotNum,:,sliceNum),3))),'',[0,1]);

% SoS on images
xtAcqShot = zeros(dim.x, dim.y, dim.dyn/ktFactor, 1, dim.loca);
for iShot = 1:(dim.dyn/ktFactor)
    for iSlice = 1:dim.loca
        xtAcqSoS(:,:,iShot,:,iSlice) = sqrt( sum( abs( sum(xtAcqSS(:,:,:,iShot,:,iSlice),3) ).^2 ,5) );
    end
end
montage_RR(squeeze(xtAcqSoS(:,:,shotNum,:,sliceNum)),'',[0,1]);

% DC SoS Image
montage_RR(squeeze( sum( xtAcqSoS(:,:,:,:,sliceNum), 3)),'',[0,50] );


%% Create Sum of Squares Data Equivalent to Body Coil Acquisition
% ktAcqSoS = reshape(ktAcq , dim.x, dim.y, ktFactor, dim.dyn/ktFactor, dim.chan, dim.loca );
% ktAcqSoS = sum(ktAcqSoS.^2 , 5);
% 
% ktAcqSoS = phase_correct( ktAcqSoS );
% xtAcqSoS = kt2xt( ktAcqSoS );
% 
% % View SoS Images
% montage_RR(abs(squeeze(sqrt(sum(xtAcqSoS(:,:,:,shotNum,:,sliceNum),3)))));


%% Calculate Shot-by-shot Coil Sensitivity Maps
for iShot = 1:(dim.dyn/ktFactor)
    for iSlice = 1:dim.loca
        for iCoil = 1:dim.chan
            csmSS(:,:,iShot,:,iSlice,iCoil) = squeeze( abs( sum( xtAcqSS( :,:,:,iShot,iCoil,iSlice ), 3) ) ) ./ xtAcqSoS(:,:,iShot,:,iSlice);
        end
    end
end

montage_RR(squeeze(csmSS(:,:,shotNum,:,sliceNum,:)),['Single-shot Coil Sens. Maps, shot' num2str(shotNum)],[0,1]); colormap('jet');

montage_RR(squeeze( sum( csmSS(:,:,:,:,sliceNum,:), 3) ),['Sum of Single-shot Coil Sens. Maps'],[0,10]); colormap('jet');

% montage_RR(squeeze( sqrt(sum( csmSS(:,:,:,:,sliceNum,:).^2, 3)) ),['Sum of Squares Single-shot Coil Sens. Maps'],[0,sqrt(10)]); colormap('jet');






