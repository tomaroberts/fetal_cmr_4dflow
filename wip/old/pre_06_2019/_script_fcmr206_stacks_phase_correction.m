sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr206';
cd(sp);

dataDir = '\data';
ktReconDir = '\ktrecon';
maskDir = '\mask';

sIDs = {'s17','s18','s20','s21'};
typeStr = {'_rlt','_dc'};
maskStr = 'lungs';

for tt = 1:2

%     for ss = 1
    for ss = 1:numel(sIDs)
        
        clear ab re im ph ph_corr cx

        ab = load_untouch_nii( [sp ktReconDir '\' sIDs{ss} typeStr{tt} '_ab.nii.gz'] );
        re = load_untouch_nii( [sp ktReconDir '\' sIDs{ss} typeStr{tt} '_re.nii.gz'] );
        im = load_untouch_nii( [sp ktReconDir '\' sIDs{ss} typeStr{tt} '_im.nii.gz'] );
        
        % initialise ph_corr nifti structure based on existing real data
        ph = re;
        ph_corr = re;
        
        cx.img = re.img + i.*im.img;
        ph.img = angle(abs(cx.img).*exp(sqrt(-1)*(angle(cx.img)+pi)));
        
        % save ph nifti (without polynomial correction)
        ph.img = single(ph.img); % convert to single to match original .nii
        save_untouch_nii(ph,[sp dataDir '\' sIDs{ss} typeStr{tt} '_ph.nii.gz']);
        
        disp(['Created ' sIDs{ss} typeStr{tt} '_ph.nii.gz ...']);
   
        mask_poly = load_untouch_nii( [sp maskDir '\' sIDs{ss} '_mask_' maskStr '.nii.gz'] );
%         mask_amnio = load_untouch_nii( [sp maskDir '\' sIDs{ss} '_mask_amniotic_fluid.nii.gz'] );
%         mask_amnio = load_untouch_nii( [sp maskDir '\' sIDs{ss} '_mask_lungs.nii.gz'] );
%         mask_heart = load_untouch_nii( [sp maskDir '\' sIDs{ss} '_mask_heart.nii.gz'] );
        
        ph_corr.img = zeros(size(ph.img));
        
        for ii = 1:size(ph.img,4)
            
            %%% NB: need to check this is actually 3D poly and not
            %%% series of slice-by-slice 2D polys
            
            ph_corr.img(:,:,:,ii) = poly3D_phase_subtract(ph.img(:,:,:,ii),mask_poly.img);
        end
        
       
%         % view
%         %-mask stacks
%         for ii = 1:size(ab.img,4); ab_m.img(:,:,:,ii) = ab.img(:,:,:,ii).*double(mask_heart.img); end
%         for ii = 1:size(ph.img,4); ph_m.img(:,:,:,ii) = ph.img(:,:,:,ii).*double(mask_heart.img); end
%         for ii = 1:size(ph_corr.img,4); ph_corr_m.img(:,:,:,ii) = ph_corr.img(:,:,:,ii).*double(mask_heart.img); end
%         
%         sliceNum = 8;
%         
%         implay_RR([ab.img(:,:,sliceNum,:), ab_m.img(:,:,sliceNum,:)]);
%         implay_RR([ph.img(:,:,sliceNum,:), ph_corr.img(:,:,sliceNum,:)],'gray',[-pi/2,pi/2]);
%         implay_RR([ph_corr.img(:,:,sliceNum,:), ph_corr_m.img(:,:,sliceNum,:), ph_m.img(:,:,sliceNum,:)],'gray',[-pi/2,pi/2]);
        
        
%         % shift phase by 1000 to make all positive
%         shiftFactor = 1000;
%         ph_corr.img = ph_corr.img + shiftFactor;
%         ph_corr.hdr.dime.glmax = min(ph_corr.img(:));
%         ph_corr.hdr.dime.glmin = max(ph_corr.img(:));
        
        % save ph_corr.nii.gz
        ph_corr.img = single(ph_corr.img); % convert to single to match original .nii
        save_untouch_nii(ph_corr,[sp dataDir '\' sIDs{ss} typeStr{tt} '_ph_corr_' maskStr '.nii.gz']);
        

        
        disp(['Created ' sIDs{ss} typeStr{tt} '_ph_corr_' maskStr '.nii.gz ...']);
    end
end



