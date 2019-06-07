%% Correcting 'grass' in fcmr214

fcmrDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr191\';
velDir = 'vel_vol_rg-nrs_phcorr_uterus';

% fcmrDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr194\';
% velDir = 'vel_vol_rg-nrs_phcorr_uterus';

% fcmrDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr214\';
% velDir = 'vel_vol-rg-nrs_phcorr_uterus';

cd(fcmrDir);

cv = load_untouch_nii('cine_vol/cine_vol.nii.gz');
v0 = load_untouch_nii([velDir '/velocity-final-0.nii.gz']);
v1 = load_untouch_nii([velDir '/velocity-final-1.nii.gz']);
v2 = load_untouch_nii([velDir '/velocity-final-2.nii.gz']);

% Resample cine if dimensions do not exactly match velocity volume
if size(cv.img,1) ~= size(v0.img,1) || size(cv.img,2) ~= size(v0.img,2) || size(cv.img,3) ~= size(v0.img,3)
    for tt = 1:size(cv.img,4)
        cv_re.img(:,:,:,tt) = imresize3( cv.img(:,:,:,tt) , size(v0.img(:,:,:,1)) );
    end
    cv.img = cv_re.img;
end

% Mask blood pool
bpThreshold = 80; % 'blood pool' threshold - arbitrary number empirically determined
BP = cv.img>bpThreshold;
% B = cv.img<prctile(cv.img(:),60); % percentile method

% View masked blood pool / bright regions
% montage_RR(cv.img(:,:,:,1),'gray',[0,150]);
% montage_RR(BP(:,:,:,1) .* cv.img(:,:,:,1),'gray',[0,150]);

% Mask background (and velocity 'rim') and in v0/v1/v2 --- currently set to -15
% BK = cv.img<0;
BK = v0.img==-15;
v0.img(BK) = 0; v1.img(BK) = 0; v2.img(BK) = 0;

% Temporal mask
% - ONLY allow voxels which never see blood pool through time
NBP = logical(~BP & ~BK);           % non-blood pool
MT = sum(NBP,4);                    % sum masks through time
MT( MT < size(cv.img,4) ) = 0;      % set any voxels which see blood pool = 0
MT = logical(MT);

% Find average velocity through time in non-blood pool regions
v0_mean_time.img = mean( v0.img , 4  ) .* MT;
v1_mean_time.img = mean( v1.img , 4  ) .* MT;
v2_mean_time.img = mean( v2.img , 4  ) .* MT;

% Find median velocity through time in non-blood pool regions
v0_median_time.img = median( v0.img , 4 ) .* MT;
v1_median_time.img = median( v1.img , 4 ) .* MT;
v2_median_time.img = median( v2.img , 4 ) .* MT;

% STDEV velocity through time in non-blood pool regions
v0_std_time.img = std( v0.img , [], 4  ) .* MT;
v1_std_time.img = std( v1.img , [], 4  ) .* MT;
v2_std_time.img = std( v2.img , [], 4  ) .* MT;


%% Simple subtraction of median/mean velocity from each timepoint

% % Subtract median velocity from each timepoint
% v0_corr.img = v0.img - v0_median_time.img;
% v1_corr.img = v1.img - v1_median_time.img;
% v2_corr.img = v2.img - v2_median_time.img;

% % Subtract average velocity from each timepoint
% v0_corr.img = v0.img - v0_mean_time.img;
% v1_corr.img = v1.img - v1_mean_time.img;
% v2_corr.img = v2.img - v2_mean_time.img;


%% Apply polynomial to DC velocity volume
polyOrder = 3;
M = logical(MT);

[xdc,ydc,zdc] = meshgrid( 1:size(cv.img,2) , 1:size(cv.img,1) , 1:size(cv.img,3) );
model0 = polyfitn( [xdc(M(:)),ydc(M(:)),zdc(M(:))], v0_median_time.img(M(:)), polyOrder );
model1 = polyfitn( [xdc(M(:)),ydc(M(:)),zdc(M(:))], v1_median_time.img(M(:)), polyOrder );
model2 = polyfitn( [xdc(M(:)),ydc(M(:)),zdc(M(:))], v2_median_time.img(M(:)), polyOrder );

% x = repmat(xdc,[1,1,nFrame]);
% y = repmat(ydc,[1,1,nFrame]);
% z = repmat(zdc,[1,1,nFrame]);
% model0 = polyfitn( [x(M(:)),y(M(:)),z(M(:))], v0_median_time.img(M(:)), polyOrder );
% model1 = polyfitn( [x(M(:)),y(M(:)),z(M(:))], v1_median_time.img(M(:)), polyOrder );
% model2 = polyfitn( [x(M(:)),y(M(:)),z(M(:))], v2_median_time.img(M(:)), polyOrder );

V0 = reshape( polyvaln( model0, [xdc(:),ydc(:),zdc(:)] ), [size(cv.img,1),size(cv.img,2),size(cv.img,3)] );
V1 = reshape( polyvaln( model1, [xdc(:),ydc(:),zdc(:)] ), [size(cv.img,1),size(cv.img,2),size(cv.img,3)] );
V2 = reshape( polyvaln( model2, [xdc(:),ydc(:),zdc(:)] ), [size(cv.img,1),size(cv.img,2),size(cv.img,3)] );

% Subtract polynomial from each timepoint
for tt = 1:size(v0.img,4)
    v0_corr.img(:,:,:,tt) = v0.img(:,:,:,tt) - V0;
    v1_corr.img(:,:,:,tt) = v1.img(:,:,:,tt) - V1;
    v2_corr.img(:,:,:,tt) = v2.img(:,:,:,tt) - V2;
end

% Apply background mask
v0_corr.img = v0_corr.img .* ~BK;
v1_corr.img = v1_corr.img .* ~BK;
v2_corr.img = v2_corr.img .* ~BK;

imtar([v0.img(:,:,16,20) V0(:,:,16) v0_corr.img(:,:,16,20); ...
       v0_median_time.img(:,:,16) v0_std_time.img(:,:,16) zeros(size(v0_std_time.img(:,:,16)))],-0.05,0.05);

% View raw and polynomial corrected velocity volumes
implay_RR([v0.img(:,:,16,:) v0_corr.img(:,:,16,:); ...
           v1.img(:,:,16,:) v1_corr.img(:,:,16,:); ...
           v2.img(:,:,16,:) v2_corr.img(:,:,16,:); ...
           v0.img(:,:,16,:)+v1.img(:,:,16,:)+v2.img(:,:,16,:) v0_corr.img(:,:,16,:)+v1_corr.img(:,:,16,:)+v2_corr.img(:,:,16,:)], ...
           'jet',[-0.2,0.2]);

%% Save background corrected .nii
v0.img = v0_corr.img; v1.img = v1_corr.img; v2.img = v2_corr.img;

% save_untouch_nii(v0,[velDir '/velocity-final-bkCorr-0.nii.gz']);
% save_untouch_nii(v1,[velDir '/velocity-final-bkCorr-1.nii.gz']);
% save_untouch_nii(v2,[velDir '/velocity-final-bkCorr-2.nii.gz']);

save_untouch_nii(v0,[velDir '/velocity-final-polyCorr-0.nii.gz']);
save_untouch_nii(v1,[velDir '/velocity-final-polyCorr-1.nii.gz']);
save_untouch_nii(v2,[velDir '/velocity-final-polyCorr-2.nii.gz']);

%% Save for MRtrix

% Flip third velocity component for MRtrix

% saveDir = [velDir '/vel_vol_4d/2flipped_bkCorr/'];
% saveName = '_2flipped-bkCorr-vel-4D-t';

saveDir = [velDir '/vel_vol_4d/2flipped_polyCorr/'];
saveName = '_2flipped-polyCorr-vel-4D-t';

mkdir(saveDir); 
v2.img = -1 .* v2.img;

for ii = 1:size(v0.img,4)
    
    %v012
    v3D = v0;
    v3D.img = cat(4,v0.img(:,:,:,ii), v1.img(:,:,:,ii), v2.img(:,:,:,ii));
    
    v3D.hdr.dime.dim(5) = 3;
    v3D.hdr.dime.pixdim(5) = 1;
    save_untouch_nii(v3D,[saveDir 'v012' saveName num2str(ii) '.nii.gz']);

end


% Default components for Paraview
saveDir = [velDir '/vel_vol_4d/default_polyCorr/'];
saveName = '_default-polyCorr-vel-4D-t';

mkdir(saveDir); 
v2.img = v2_corr.img;

for ii = 1:size(v0.img,4)
    
    %v012
    v3D = v0;
    v3D.img = cat(4,v0.img(:,:,:,ii), v1.img(:,:,:,ii), v2.img(:,:,:,ii));
    
    v3D.hdr.dime.dim(5) = 3;
    v3D.hdr.dime.pixdim(5) = 1;
    save_untouch_nii(v3D,[saveDir 'v012' saveName num2str(ii) '.nii.gz']);

end