%% 2018_11_01 --- PC-bSSFP --- FlowPhantom_6 --- spherical flask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP';
cd(sp);
cd affine_registration

%% Step 1)
%% Files reconstructed on beastie01 using mrecon_pc_recon().m
%- outputs _CPLX / _MAG / _PH.nii.gz files

% % rawFileName === array of .raw filenames.
% % senseRefFileName = 'fl_01112018_1959030_1000_11_pih1senserefscanclearV4.lab';
% % coilSurveyFileName = 'fl_01112018_1958375_1000_8_pih1coilsurveyscanV4.lab';
% % senseflag = 'on';
% % reconArray = 'classical';
% % 
% % for ii = 1:numel(rawFileName)
% %     mrecon_pc_recon( rawFileName{ii},...
% %         'reconArray', reconArray,...
% %         'senseref', senseRefFileName,...
% %         'coilsurvey', coilSurveyFileName,...
% %         'senseflag', senseflag );
% %     
% %     disp(['Reconned file ' num2str(ii) '...']);
% % end


%% Step 2)
%% Transform stacks into same space using irtk transformation tool
% - target stack = tra_plus_MAG (could be any)
% - can transform PHASE data to tra_plus_MAG too
% - code on gpubeastie: 
%   transformation sag_plus_MAG.nii.gz sag_plus_MAG_transformed.nii.gz -target tra_plus_MAG.nii.gz -linear -matchInputType
% - repeat for all stacks
% - also transformed static_water_mask.nii.gz (which was drawn in MITK)


%% Step 3)
%% Affine transformation using irtk
% - - target stack = tra_plus_MAG_transformed (could be any)
% - code on gpubeastie: 
%   areg tra_plus_MAG_transformed.nii.gz sag_plus_MAG_transformed.nii.gz -dofout areg_tp_sp.dof
%   transformation sag_plus_MAG_transformed.nii.gz sag_plus_MAG_affine.nii.gz -dofin areg_tp_sp.dof -target tra_plus_MAG_transformed.nii.gz -linear
% - repeat for all stacks


%% Step 3)
%% Crop common area

%% LINEAR TRANSFORMATION ONLY
% nii = load_untouch_nii('tra_plus_MAG_transformed.nii.gz'); MAG.tra.plus = nii.img; 
% nii = load_untouch_nii('sag_plus_MAG_transformed.nii.gz'); MAG.sag.plus = nii.img; 
% nii = load_untouch_nii('cor_plus_MAG_transformed.nii.gz'); MAG.cor.plus = nii.img; 
% nii = load_untouch_nii('tra_minus_MAG_transformed.nii.gz'); MAG.tra.minus = nii.img; 
% nii = load_untouch_nii('sag_minus_MAG_transformed.nii.gz'); MAG.sag.minus = nii.img; 
% nii = load_untouch_nii('cor_minus_MAG_transformed.nii.gz'); MAG.cor.minus = nii.img;
% 
% nii = load_untouch_nii('tra_plus_PH_transformed.nii.gz'); PH.tra.plus = nii.img; 
% nii = load_untouch_nii('sag_plus_PH_transformed.nii.gz'); PH.sag.plus = nii.img; 
% nii = load_untouch_nii('cor_plus_PH_transformed.nii.gz'); PH.cor.plus = nii.img; 
% nii = load_untouch_nii('tra_minus_PH_transformed.nii.gz'); PH.tra.minus = nii.img; 
% nii = load_untouch_nii('sag_minus_PH_transformed.nii.gz'); PH.sag.minus = nii.img; 
% nii = load_untouch_nii('cor_minus_PH_transformed.nii.gz'); PH.cor.minus = nii.img;

% % sagittal stack shifted by one compared to others, so view post-shift:
% implay_RR([MAG.tra.plus(:,:,1:end-1), MAG.sag.plus(:,:,2:end), MAG.cor.plus(:,:,1:end-1);...
%            MAG.tra.minus(:,:,1:end-1), MAG.sag.minus(:,:,2:end), MAG.cor.minus(:,:,1:end-1)]);
       
% % manually determined x/y crop values by looking at full stacks
% xcrop = 109:180;
% ycrop = 134:197;
% zcrop_s = 2:size(MAG.sag.plus,3);
% zcrop_t_c = 1:size(MAG.tra.plus,3)-1;
%        
% % apply crop to stacks
% MAG.tra.plus = MAG.tra.plus(xcrop,ycrop,zcrop_t_c);   PH.tra.plus = PH.tra.plus(xcrop,ycrop,zcrop_t_c);
% MAG.sag.plus = MAG.sag.plus(xcrop,ycrop,zcrop_s);   PH.sag.plus = PH.sag.plus(xcrop,ycrop,zcrop_s);
% MAG.cor.plus = MAG.cor.plus(xcrop,ycrop,zcrop_t_c);   PH.cor.plus = PH.cor.plus(xcrop,ycrop,zcrop_t_c);  
% MAG.tra.minus = MAG.tra.minus(xcrop,ycrop,zcrop_t_c); PH.tra.minus = PH.tra.minus(xcrop,ycrop,zcrop_t_c);
% MAG.sag.minus = MAG.sag.minus(xcrop,ycrop,zcrop_s); PH.sag.minus = PH.sag.minus(xcrop,ycrop,zcrop_s);
% MAG.cor.minus = MAG.cor.minus(xcrop,ycrop,zcrop_t_c); PH.cor.minus = PH.cor.minus(xcrop,ycrop,zcrop_t_c);
% MASK = single(MASK(xcrop,ycrop,zcrop_t_c));
% 
% clear xcrop ycrop zcrop_s zcrop_t_c



%% AFFINE APPROACH
%-cartesian stacks
nii = load_untouch_nii('tra_plus_MAG_affine.nii.gz'); MAG.tra.plus = nii.img; 
nii = load_untouch_nii('sag_plus_MAG_affine.nii.gz'); MAG.sag.plus = nii.img; 
nii = load_untouch_nii('cor_plus_MAG_affine.nii.gz'); MAG.cor.plus = nii.img; 
nii = load_untouch_nii('tra_minus_MAG_affine.nii.gz'); MAG.tra.minus = nii.img; 
nii = load_untouch_nii('sag_minus_MAG_affine.nii.gz'); MAG.sag.minus = nii.img; 
nii = load_untouch_nii('cor_minus_MAG_affine.nii.gz'); MAG.cor.minus = nii.img;

nii = load_untouch_nii('tra_plus_PH_affine.nii.gz'); PH.tra.plus = nii.img; Srow.tra.plus = getSrow(load_untouch_nii('tra_plus_PH.nii.gz')); %<-get Srow from original .nii
nii = load_untouch_nii('sag_plus_PH_affine.nii.gz'); PH.sag.plus = nii.img; Srow.sag.plus = getSrow(load_untouch_nii('sag_plus_PH.nii.gz'));
nii = load_untouch_nii('cor_plus_PH_affine.nii.gz'); PH.cor.plus = nii.img; Srow.cor.plus = getSrow(load_untouch_nii('cor_plus_PH.nii.gz'));
nii = load_untouch_nii('tra_minus_PH_affine.nii.gz'); PH.tra.minus = nii.img; Srow.tra.minus = getSrow(load_untouch_nii('tra_minus_PH.nii.gz'));
nii = load_untouch_nii('sag_minus_PH_affine.nii.gz'); PH.sag.minus = nii.img; Srow.sag.minus = getSrow(load_untouch_nii('sag_minus_PH.nii.gz'));
nii = load_untouch_nii('cor_minus_PH_affine.nii.gz'); PH.cor.minus = nii.img; Srow.cor.minus = getSrow(load_untouch_nii('cor_minus_PH.nii.gz'));

nii = load_untouch_nii('static_water_mask_transformed.nii.gz'); MASK = nii.img; %drawn on MAG_tra_minus

clear nii

% manually determined x/y crop values by looking at full stacks
xcrop = 109:180;
ycrop = 134:197;
       
% apply crop to stacks
MAG.tra.plus = MAG.tra.plus(xcrop,ycrop,:);   PH.tra.plus = PH.tra.plus(xcrop,ycrop,:);
MAG.sag.plus = MAG.sag.plus(xcrop,ycrop,:);   PH.sag.plus = PH.sag.plus(xcrop,ycrop,:);
MAG.cor.plus = MAG.cor.plus(xcrop,ycrop,:);   PH.cor.plus = PH.cor.plus(xcrop,ycrop,:);  
MAG.tra.minus = MAG.tra.minus(xcrop,ycrop,:); PH.tra.minus = PH.tra.minus(xcrop,ycrop,:);
MAG.sag.minus = MAG.sag.minus(xcrop,ycrop,:); PH.sag.minus = PH.sag.minus(xcrop,ycrop,:);
MAG.cor.minus = MAG.cor.minus(xcrop,ycrop,:); PH.cor.minus = PH.cor.minus(xcrop,ycrop,:);
MASK = single(MASK(xcrop,ycrop,:));

% view stacks
% implay_RR([MAG.tra.plus, MAG.sag.plus, MAG.cor.plus;...
%            MAG.tra.minus, MAG.sag.minus, MAG.cor.minus]);
            
% implay_RR([PH.tra.plus, PH.sag.plus, PH.cor.plus;...
%            PH.tra.minus, PH.sag.minus, PH.cor.minus]);

% view bipolar deltaPhase images
% implay_RR([PH.tra.plus-PH.tra.minus,
%            PH.sag.plus-PH.sag.minus,
%            PH.cor.plus-PH.cor.minus],'jet',[-pi,pi]);

% % view normalised, subtracted MAG images
% maxtra = max(MAG.tra.plus(:)); MAG.tra.plus = MAG.tra.plus./maxtra;
% maxsag = max(MAG.sag.plus(:)); MAG.sag.plus = MAG.sag.plus./maxsag;
% maxcor = max(MAG.cor.plus(:)); MAG.cor.plus = MAG.cor.plus./maxcor;
% implay_RR([MAG.tra.plus, MAG.sag.plus, MAG.cor.plus]);


%-oblique stacks
nii = load_untouch_nii('obl_tra_ap-45_MAG_affine.nii.gz'); MAG.obl.tra = nii.img; 
nii = load_untouch_nii('obl_sag_lr45_MAG_affine.nii.gz'); MAG.obl.sag = nii.img; 
nii = load_untouch_nii('obl_cor_ap45_MAG_affine.nii.gz'); MAG.obl.cor = nii.img;

nii = load_untouch_nii('obl_tra_ap-45_PH_affine.nii.gz'); PH.obl.tra = nii.img; Srow.obl.tra = getSrow(load_untouch_nii('obl_tra_ap-45_PH.nii.gz'));
nii = load_untouch_nii('obl_sag_lr45_PH_affine.nii.gz'); PH.obl.sag = nii.img; Srow.obl.sag = getSrow(load_untouch_nii('obl_sag_lr45_PH.nii.gz'));
nii = load_untouch_nii('obl_cor_ap45_PH_affine.nii.gz'); PH.obl.cor = nii.img; Srow.obl.cor = getSrow(load_untouch_nii('obl_cor_ap45_PH.nii.gz'));

% implay_RR([MAG.obl.tra, MAG.obl.sag, MAG.obl.cor]);

% apply crop to stacks
MAG.obl.tra = MAG.obl.tra(xcrop,ycrop,:);   PH.obl.tra = PH.obl.tra(xcrop,ycrop,:);
MAG.obl.sag = MAG.obl.sag(xcrop,ycrop,:);   PH.obl.sag = PH.obl.sag(xcrop,ycrop,:);
MAG.obl.cor = MAG.obl.cor(xcrop,ycrop,:);   PH.obl.cor = PH.obl.cor(xcrop,ycrop,:);  

% tidy workspace
clear xcrop ycrop zcrop_s zcrop_t_c nii


%% Step 4)
%% Correct phase rolls in stacks

%-cartesian stacks
Y0mat(:,:,:,1) = PH.tra.plus;
Y0mat(:,:,:,2) = PH.sag.plus;
Y0mat(:,:,:,3) = PH.cor.plus;
Y0mat(:,:,:,4) = PH.tra.minus;
Y0mat(:,:,:,5) = PH.sag.minus;
Y0mat(:,:,:,6) = PH.cor.minus;

%-oblique stacks
Y0mat(:,:,:,7) = PH.obl.tra;
Y0mat(:,:,:,8) = PH.obl.sag;
Y0mat(:,:,:,9) = PH.obl.cor;

for ii = 1:size(Y0mat,4)
    Y0 = Y0mat(:,:,:,ii);

    M = MASK;
    [nRow,nCol,nSl,nFrame] = size( Y0 );

    % EstimateConstant Background Phase
    P0 = median(Y0(find(M)));

    % Fit Polynomial to Background Phase
    [xdc,ydc,zdc] = meshgrid(1:nCol,1:nRow,1:nSl);
    % [xdc,ydc,zdc] = meshgrid(voxSize(2)*(1:nCol),voxSize(1)*(1:nRow),voxSize(3)*(1:nSl));
    x = repmat(xdc,[1,1,nFrame]);
    y = repmat(ydc,[1,1,nFrame]);
    z = repmat(zdc,[1,1,nFrame]);
    polyOrder = 3;
    % model = polyfitn( [x(find(M)),y(find(M)),z(find(M))], angle(Y0(find(M))*exp(-(1i*P0))), polyOrder );
    model = polyfitn( [x(find(M)),y(find(M)),z(find(M))], Y0(find(M)), polyOrder );

    P1 = reshape( polyvaln( model, [xdc(:),ydc(:),zdc(:)] ), [nRow,nCol,nSl] );

    Y(:,:,:,ii) = Y0-P1;
end

% view phase corrected stacks
implay_RR([Y(:,:,:,1), Y(:,:,:,2), Y(:,:,:,3);...
           Y(:,:,:,4), Y(:,:,:,5), Y(:,:,:,6);...
           Y(:,:,:,7), Y(:,:,:,8), Y(:,:,:,9)]);
       
% view bipolar corrected deltaPhase images
implay_RR([Y(:,:,:,1)-Y(:,:,:,4), Y(:,:,:,2)-Y(:,:,:,5), Y(:,:,:,3)-Y(:,:,:,6)],'jet',[-pi,pi]);


%% Step 5)
%% Make bipolar velocity maps

Ms = -6.65;
VELbip.venc = getVenc_bipolar(Ms);

Ms = -6.65 .* (1e-3).^2; %<-bSSFP
% Ms = 5.88  .* (1e-3).^2; %<- QFlow with VENC=100cm/s
gamma = 2 .* pi .* 42577;

% venc = pi ./ ( 2 .* gamma .* Ms ); %<- m/s
% VELbip.venc = venc .* 1e2; %<- cm/s

VELbip.tra = ( Y(:,:,:,1)-Y(:,:,:,4) ) ./ ( 2 .* gamma .* Ms ); VELbip.tra = VELbip.tra .* 1e2; %<- cm/s
VELbip.sag = ( Y(:,:,:,2)-Y(:,:,:,5) ) ./ ( 2 .* gamma .* Ms ); VELbip.sag = VELbip.sag .* 1e2; %<- cm/s
VELbip.cor = ( Y(:,:,:,3)-Y(:,:,:,6) ) ./ ( 2 .* gamma .* Ms ); VELbip.cor = VELbip.cor .* 1e2; %<- cm/s

implay_RR([VELbip.tra, VELbip.sag, VELbip.cor],'jet',[-abs(VELbip.venc),abs(VELbip.venc)]);


%% Step 6)
%% Convert referenceless phase maps to velocity maps

clear P Pmat V VEL

Mm =  6.21 .* (1e-3).^2;  %ms^2.mT/m
Mp =  0;
Ms = -6.65 .* (1e-3).^2;  %ms^2.mT/m
gamma = 1.5 .* 2 .* pi .* 42577; %Hz/mT --- 1.5 because 1.5T
% gamma = 2 .* pi .* 42577; %Hz/mT

%-------LR--AP-FH
Mtra = [Mm, 0, Ms];
Msag = [Ms, 0, Mm];
% Msag = [Ms, Mm, 0]; %from Srow
Mcor = [Mm, Ms, 0];

%flipped bSSFP scans
Mtra_ = [Mm, 0, -Ms];
Msag_ = [-Ms, 0, Mm];
% Msag_ = [-Ms, Mm, 0]; %from Srow
Mcor_ = [Mm, -Ms, 0];

%oblique bSSFP scans
Mtra_ap45_ = [Mm*cos(d2r(45)) 0 -Ms*cos(d2r(45))];  %rotate around AP, PE dir stays same
% Mtra_ap45_ = [Mm*cos(d2r(45))+Ms*cos(d2r(45)) 0 -Ms*cos(d2r(45))+Mm*cos(d2r(45))];
Msag_lr45 =  [Ms Mm*cos(d2r(45)) Mm*cos(d2r(45))];   %rotate around LR, slice dir stays same
Mcor_ap45 =  [Mm*cos(d2r(45)), Ms, Mm*cos(d2r(45))]; %rotate around AP, slice dir stays same

% vectorise phase images
P.tra = Y(:,:,:,1); P.tra = P.tra(:);
P.sag = Y(:,:,:,2); P.sag = P.sag(:);
P.cor = Y(:,:,:,3); P.cor = P.cor(:);
P.tra_ = Y(:,:,:,4); P.tra_ = P.tra_(:);
P.sag_ = Y(:,:,:,5); P.sag_ = P.sag_(:);
P.cor_ = Y(:,:,:,6); P.cor_ = P.cor_(:);

P.tra_ap45_ = Y(:,:,:,7); P.tra_ap45_ = P.tra_ap45_(:);
P.sag_lr45 = Y(:,:,:,8); P.sag_lr45 = P.sag_lr45(:);
P.cor_ap45 = Y(:,:,:,9); P.cor_ap45 = P.cor_ap45(:);

% solve
% M6 = gamma .* [Mtra; Msag; Mcor; Mtra_; Msag_; Mcor_];
% Pmat = [P.tra, P.sag, P.cor, P.tra_, P.sag_, P.cor_]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M6\Pmat(:,ii);
% end

% M6obl = gamma .* [Mtra; Msag; Mcor; -Mtra_ap45_; Msag_lr45; Mcor_ap45];
% Pmat = [P.tra, P.sag, P.cor, P.tra_ap45_, P.sag_lr45, P.cor_ap45]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M6obl\Pmat(:,ii);
% end

% M5 = gamma .* [Mtra; Msag; Mcor; Mtra_; Msag_];
% Pmat = [P.tra, P.sag, P.cor, P.tra_, P.sag_]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M5\Pmat(:,ii);
% end

% M5obl = gamma .* [Mtra; Msag; Mcor; -Mtra_ap45_; Mcor_ap45];
% Pmat = [P.tra, P.sag, P.cor, P.tra_ap45_, P.cor_ap45]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M5\Pmat(:,ii);
% end

% M4 = gamma .* [Mtra; Msag; Mtra_; Msag_];
% Pmat = [P.tra, P.sag, P.cor, P.tra_]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M4\Pmat(:,ii);
% end

% M4 = gamma .* [Mtra; Msag; Mtra_; Mcor_];
% Pmat = [P.tra, P.sag, P.cor, P.cor_]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M4\Pmat(:,ii);
% end

% M4obl = gamma .* [Mtra; Mcor; Msag_lr45; Mcor_ap45];
% Pmat = [P.tra, P.cor, P.sag_lr45, P.cor_ap45]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M4obl\Pmat(:,ii);
% end

% M4obl = gamma .* [Mtra; -Mtra_ap45_; Msag_lr45; Mcor_ap45]; %% Best yet...
M4obl = gamma .* [Mtra; -Mtra_ap45_; Msag_lr45; Mcor_ap45];
Pmat = [P.tra, P.tra_ap45_, P.sag_lr45, P.cor_ap45]';
for ii = 1:length(Pmat)
    V(:,ii) = M4obl\Pmat(:,ii);
end

% M3 = gamma .* [Mtra; Msag; Mcor];
% Pmat = [P.tra, P.sag, P.cor]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M3\Pmat(:,ii);
% end

% M3obl = gamma .* [Mtra; Msag_lr45; Mcor_ap45];
% Pmat = [P.tra, P.sag_lr45, P.cor_ap45]';
% for ii = 1:length(Pmat)
%     V(:,ii) = M3obl\Pmat(:,ii);
% end

% un-vectorise and convert to velocity stacks
subV=ind2subV(size(PH.tra.plus),1:numel(P.tra));

for ii = 1:length(subV)
    i = subV(ii,1);
    j = subV(ii,2);
    k = subV(ii,3);
    VEL.vx(i,j,k) = V(1,ii) .* 1e2; %cm/s
    VEL.vy(i,j,k) = V(2,ii) .* 1e2;
    VEL.vz(i,j,k) = V(3,ii) .* 1e2;
end

disp('VEL done');

% view

% view referenceless velocity maps
% implay_RR([VEL.vz, VEL.vx, VEL.vy],'jet',[-50,50]);

% view comparison with bipolar maps
%%% This seems to work for Cartesian! ie: M6
% implay_RR(([VELbip.sag, VELbip.cor, VELbip.tra;...
%            VEL.vx, VEL.vy, VEL.vz]),'jet',[-50,50]);

%%% This seems to work for my current combination of oblique! ie: M4obl
% implay_RR(([VELbip.sag, VELbip.cor, VELbip.tra;...
%            VEL.vx, VEL.vy, VEL.vz]),'jet',[-50,50]);
       
% % view magnitude of 3-vectors
% implay_RR(([sqrt(VELbip.sag.^2 + VELbip.cor.^2 + VELbip.tra.^2);...
%            sqrt(VEL.vx.^2 + VEL.vy.^2 + VEL.vz.^2)]),'jet',[-50,50]);       
       

% QFLOW vs. PC-bSSFP Bipolar vs. Referenceless PC-bSSFP
load('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\QFlow\affine_registration\FlowPhantom_QFlow_4D.mat');
implay_RR(([QFLOW.VEL.sag, QFLOW.VEL.cor, QFLOW.VEL.tra;...
           -VELbip.sag/2, VELbip.cor/2, VELbip.tra/2;...
           ...-VEL.vx, VEL.vy, VEL.vz ...
           ]),'jet',[-50,50]);

       
%% save data
PCbSSFP.MAG  = MAG;
PCbSSFP.PH   = Y; %<- phase roll corrected
PCbSSFP.VEL  = VEL;
PCbSSFP.VELbip  = VELbip;
PCbSSFP.MASK = MASK;
save('FlowPhanton_PCbSSFP_4D.mat','PCbSSFP');

       
%% Make GIFs
frms = size(VEL.vx,3);
figure();

for tt = 1:frms
    imagesc([VELbip.tra(:,:,tt), VELbip.sag(:,:,tt), -VELbip.cor(:,:,tt);...
             VEL.vz(:,:,tt), VEL.vx(:,:,tt), VEL.vy(:,:,tt)],[-50,50]);
    axis image; colormap(jet);
    truesize([600 600]);
    set(gcf,'defaultfigurecolor',[1 1 1]); % remove grey background?
    set(gca,'visible','off');
    f = getframe(gcf);
    imgif{tt} = frame2im(f);
    clf;
end
% write_gif_lj(imgif,[pwd '\PC-bSSFP_bipolar_v_referenceless.gif'],'DelayTime',0.15);
write_gif_lj(imgif,[pwd '\PC-bSSFP_bipolar_v_referenceless_4stacks_tra_3obl.gif'],'DelayTime',0.15);
close all




