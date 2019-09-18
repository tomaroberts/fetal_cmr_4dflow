cd('E:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\4D_Flow_Paper\fcmr254_tom_recon\ktrecon_s20_slice8_adaptive');

xt   = niftiread('s20_rlt_ab.nii.gz');
info = niftiinfo('s20_rlt_ab.nii.gz');

sliceNum = 8;

imtar(xt(:,:,sliceNum,1));
implay_RR(squeeze(xt(:,:,sliceNum,:)));

dimT = 4;

xt2xf = @( xt ) fftshift( fft( xt, [], dimT ), dimT );

xf = xt2xf(xt);

implay_RR(abs(squeeze(xf(:,:,sliceNum,:))));

xfLine  = 120;
xfPrior = squeeze(abs(xf(xfLine,:,sliceNum,:)));
xfPrior(:,49) = 0; %TODO: automatically find DC
imtar(xfPrior,0,100); title(['x-f, slice ' num2str(sliceNum)]);
imtar(xfPrior(heartBounds,:),0,500); title('Masked Region');

heartBounds = 52:70;

xfPriorMean = mean(xfPrior,1);
figure; plot(xfPriorMean);
title(['mean(x)-f, slice ' num2str(sliceNum)]);

xfPriorMean = mean(xfPrior(heartBounds,:),1);
figure; plot(xfPriorMean);
title('mean(x)-f in Masked Region');
