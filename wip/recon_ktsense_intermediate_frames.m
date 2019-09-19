%% Create Intermediate frames from kt-SENSE data

%% Admin
reconDir = 'E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr254_tom_recon';
ktDir = '\ktrecon';
scanNum = 20;
ktFactor = 8;
frameFactor = 3;
sliceNum = 8;


%% Load data
cd([reconDir ktDir]);
load(['s' num2str(scanNum) '_csm.mat']);
load(['s' num2str(scanNum) '_kspace.mat'],'ktAcq');


%% Data Dimensions
dimX = 1;  nX = size( ktAcq, dimX );    % frequency-encode
dimY = 2;  nY = size( ktAcq, dimY );    % phase-encode
dimT = 3;  nT = size( ktAcq, dimT );    % time
dimC = 4;  nC = size( ktAcq, dimC );    % channel
dimF = dimT;  nF = nT;


%% Anon Fns
kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, dimX ), dimY ) );
phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [+1;-1], size(k,1)/2, 1 ) ) ) );
inv_phase_correct = @( k ) abs(k) .* exp( sqrt(-1) * ( angle(k) + bsxfun( @times, pi/2 * ones( size(k) ), repmat( [-1;+1], size(k,1)/2, 1 ) ) ) );


%% kt Sampling Patterns
ktSmp = single( sum( sum( ktAcq, dimC ), dimX ) ~= 0 ); %spatiotemporal optimum as in Tsao paper.

% Intermediate frame kt pattern
ktSmpIntMed = sum( reshape( ktSmp, size(ktSmp,1), size(ktSmp,2), frameFactor, size(ktSmp,3)/frameFactor, size(ktSmp,4), size(ktSmp,5) ), 3);
ktSmpIntMed = permute( ktSmpIntMed, [1,2,4,5,6,3] ); % remove 3rd dimension introduced by reshape


%% Reduce Data to Single-slice, in keeping with Josh
ktAcqSlice = ktAcq(:,:,:,:,sliceNum);

% Fourier Transform
ktAcqSlice = phase_correct(ktAcqSlice); xtAcqSlice = kt2xt(ktAcqSlice);

% View Intermediate k-space frames
figure;
for ii = 1:frameFactor
    subplot(1,3,ii);
    frameLines = ii*frameFactor-frameFactor+1:ii*frameFactor;
    imagesc(abs(sum(ktAcqSlice(:,:,frameLines,1),3)),[0,100]);
    colormap(gray), colorbar('off'), axis image;
    title(['Intermediate Frame ' num2str(ii)]);
end

% View aliased Intermediate frames
figure;
for ii = 1:frameFactor
    subplot(1,3,ii);
    frameLines = ii*frameFactor-frameFactor+1:ii*frameFactor;
    imagesc(abs(sum(xtAcqSlice(:,:,frameLines,1),3)),[0,0.1]);
    colormap(gray), colorbar('off'), axis image;
    title(['Intermediate Frame ' num2str(ii)]);
end


%% Build prior images for Tikhonov Regularization
ktPri = sum(ktAcqSlice(:,:,1:8,:),3); % Fully-sampled frame
xtPri = kt2xt(ktPri);

% imtar(abs(xtPri(:,:,:,1)),0,0.5);

% SoS image
xtPriSoS = sqrt( sum( abs(xtPri).^2 ,4 ) );
imtar(abs(xtPriSoS),0,5);

% Conjugate image - gives NaNs, so not good.
load(['s' num2str(scanNum) '_csm.mat'],'csm');
csm = csm(:,:,:,:,sliceNum);

xtPriConj = sum(xtPri.*conj(csm),4)./sum(csm.*conj(csm),4);
imtar(abs(xtPriConj),0,0.5);

% SoSConj image
% - NaNs in Conj image problematic for SENSE recon, therefore take the
% noise from SoS image and apply to background of Conj image
% - not sure this is legit...
SoSbkrd = xtPriSoS <= 0.5;
Conjfgrd = ~SoSbkrd;

xtPriSoSConj = zeros(size(xtPriSoS));
xtPriSoSConj(SoSbkrd) = xtPriSoS(SoSbkrd);
xtPriSoSConj(Conjfgrd) = sqrt( xtPriConj(Conjfgrd) );
imtar(abs(xtPriSoSConj),0,1);


%% SENSE reconstruction with conventional sampling pattern
% Normal R=2 recon
% Reduced FOV version (based on R)

R = 2;

kt1 = sum(ktAcqSlice(:,:,1:R:8,:),3);
xt1 = kt2xt(kt1);

% % estimate new noise covariance
% noiseData = reshape( xt1([1,nX],:,:,:), [], nC ) * sqrt( nX * nY * nT );
% noiseCov = cov( noiseData );
psi = psiSens;

load(['s' num2str(scanNum) '_csm.mat'],'csm');
csm = csm(:,:,:,:,sliceNum);

xt1 = permute(xt1,[2,1,3,4]);
csm = permute(csm,[2,1,3,4]);

% Trim FOV as if acquired on scanner
xt1 = xt1( 1:size(xt1,1)/R ,:,:,:); 

[Nx,Ny,Nz,Nc] = size(xt1);
imReconstructed = zeros(Nx*R,Ny);

tol = 1e-9;
inv_psi = inv(psi);
% psisqinv=inv(sqrtm(psi)); % Lucilio: rather than doing S'*psi*S etc., suggested a better way of doing the calculation

for x = 1:Nx
        
    x_idx = x:Nx:Nx*R;
    
    for y = 1:Ny
        
        S = transpose(reshape(csm(x_idx,y,1,:),R,[]));
        
        unmix = pinv( S'*inv_psi*S , tol ) * S' * inv_psi;
        
%         S = S * psisqinv; % Lucilio
%         unmix = pinv(S,tol);
        
        imReconstructed(x_idx,y) = unmix*reshape(xt1(x,y,1,:),[],1);
    end
    
end

% imReconstructed(abs(imReconstructed)<1e-8)=NaN; % turns bodged voxels into NaNs
imtar(abs(imReconstructed),0,0.5);



%% SENSE reconstruction using intermediate frames
% R=3 intermediate frame recon
% Half FOV version

R = 3;

kt1 = sum(ktAcqSlice(:,:,1:3,:),3);
xt1 = kt2xt(kt1);

% % estimate new noise covariance
% noiseData = reshape( xt1([1,nX],:,:,:), [], nC ) * sqrt( nX * nY * nT );
% noiseCov = cov( noiseData );
psi = psiSens;

load(['s' num2str(scanNum) '_csm.mat'],'csm');
csm = csm(:,:,:,:,sliceNum);

xt1 = permute(xt1,[2,1,3,4]);
csm = permute(csm,[2,1,3,4]);

% Trim FOV as if acquired on scanner
xt1 = xt1(1:50,:,:,:);

% Trim prior FOV
xtPri1 = xtPriSoS;
% xtPri1 = xtPriSoSConj;
xtPri1 = permute(xtPri1,[2,1,3,4]);
xtPri1 = xtPri1(1:150,:,:,:);

% Trim csm so divisible by frameFactor
csm = csm(1:150,:,:,:);

[Nx,Ny,Nz,Nc] = size(xt1);
imReconstructed = zeros(Nx,Ny);

inv_psi = inv(psi);
tol = 1e-9;

% specify regularization method
regul = 'tik';

switch regul
    
    case 'off'

        for x = 1:Nx

            x_idx = x:Nx:Nx*R;

            for y = 1:Ny

                S = transpose(reshape(csm(x_idx,y,1,:),R,[]));
                unmix = pinv( S'*inv_psi*S , tol ) * S' * inv_psi;

        %         G = sqrt( pinv( S'*inv_psi*S ) * S'*inv_psi*S );

        %         S = S * psisqinv; % Lucilio
        %         unmix = pinv(S,tol);

                imReconstructed(x_idx,y) = unmix*reshape(xt1(x,y,1,:),[],1);
        %         gFactor(x_idx,y) = G;
            end
        end
    
    case 'tik'
        
%         lambda = [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9];
        lambda = logspace(-3,1,20);
%         lambda = 1e-1;
        
        A = eye(R);
        
        for ll = 1:numel(lambda)
            for x = 1:Nx

                x_idx = x:Nx:Nx*R;

                for y = 1:Ny

                    S = transpose(reshape(csm(x_idx,y,1,:),R,[]));
                    unmix = pinv( S'*inv_psi*S + lambda(ll).*A'*A  , tol ) * S' * inv_psi;

                    diffImPri = reshape(xt1(x,y,1,:),[],1) - S*xtPri1(x_idx,y);
                    imReconstructed(x_idx,y,ll) = xtPri1(x_idx,y) + unmix * diffImPri;               

                end
            end
            
        disp(['Finished recon no.: ' num2str(ll)]);
            
        end
        
end

montage_RR(abs(imReconstructed),'',[0,1]);
montage_RR(abs(imReconstructed)-abs(xtPri1),'diff',[-1,1]); colormap('jet');

% imtar(abs(imReconstructed),0,0.5);


