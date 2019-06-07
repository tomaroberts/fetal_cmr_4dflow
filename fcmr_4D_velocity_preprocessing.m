%% 4D Cardiac Flow - fcmr pre-processing


%% paths

% ingenia-raw
rawDir = 'Z:\';

% Josh's store:
% fcmrJoshDir = 'F:\2018_11_16_Josh_Backup\fcmr_cine_3d';
% cd(fcmrJoshDir);

% fNums = [189 194 197 201 202 213 214 230 254 255 257];
fNums = 254; % pre-process single case at a time...

% Tom local folder
fcmrTomDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data';


%% Manually draw uterus ROIs for each stack for polynomial phase correction
%- save in ../mask/ as s*_mask_uterus.nii.gz


%% Create _ph.nii.gz and _ph_corr.nii.gz

%%% TO DO?: Do I need to update the nifti headers when I make the
%%% _ph.nii.gz files? Currently headers based on _re.nii.gz images


cd(fcmrTomDir);

for ii = 1:numel(fNums)
    
    disp(['Generating phase stacks for fcmr' num2str(fNums(ii)) ' ... ' ]);
    
    cd([fcmrTomDir '\fcmr' num2str(fNums(ii)) ]);
    
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
        % FIXME: polynomials currently saving as 5D .nii - should be 3D
        poly_nii = re;
        polyNom = exp( -( 1i*(P0+P1) ) );
        
        % TODO: understand exactly what the +pi is doing. Cycling the phase
        % somehow.
        poly_nii.img = angle(abs(polyNom).*exp(sqrt(-1)*(angle(polyNom)+pi)));
        
        save_untouch_nii(poly_nii,[sIDs(ss).name(1:3) '_rlt_ph_polynomial_uterus.nii.gz']);        
        movefile([sIDs(ss).name(1:3) '_rlt_ph_polynomial_uterus.nii.gz'] , ['../data/' sIDs(ss).name(1:3) '_rlt_ph_polynomial_uterus.nii.gz']); 
        
        disp(['Created ' sIDs(ss).name(1:3) '_rlt_ph_corr_uterus.nii.gz']);
        
        
    end
    
    disp(['Completed making phase stacks for fcmr' num2str(fNums(ii)) ' ... ' ]);
    
    cd([fcmrTomDir '\fcmr' num2str(fNums(ii)) ]);
    
end


%% get scan dates

% get scan dates from Josh's log files
for ii = 1:numel(fNums)
    cd([fcmrTomDir '/fcmr' num2str(fNums(ii)) '/ktrecon/']);
    logFile = dir('log*.txt');
    
    fid = fopen(logFile(1).name);
    tline = fgetl(fid);
    while ischar(tline)
        if strfind(tline, 'pnraw01-ingenia');
            C = tline;
            break;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
    D{ii} = C(34+16:34+25); %get date
    N{ii} = C(34+27:34+35); %get id
    
    clear tline fid
end

%% Create new goalc.txt files which include ORIENT object

% 1) get raw files from pnraw

for ii = 1:numel(fNums)
    
    disp(['Fetching .raw files for fcmr' num2str(fNums(ii)) ' ... ']);
    
    cd(rawDir);
    cd(D{ii});
    cd(N{ii});
    
    rtRawScanNames = dir('*kt8i0bffe*.raw');
    rtLabScanNames = dir('*kt8i0bffe*.lab');
    
    rawDirTom = [fcmrTomDir '\fcmr' num2str(fNums(ii)) '\raw'];
    mkdir(rawDirTom);
    for nn = 1:numel(rtRawScanNames)
        copyfile(rtRawScanNames(nn).name , rawDirTom );
        copyfile(rtLabScanNames(nn).name , rawDirTom );
        disp(['Copied .raw/.lab files number ... ' num2str(nn)]);
    end
    
    cd(rawDirTom);
    
    disp(['Copied .raw/.lab files for fcmr' num2str(fNums(ii)) ' ... ']);
end


%% 2) Manually copy the ../raw folder across to beastie02
%    (if entire fcmrXXX folder is not already there)
warning('REMEMBER TO SEND ../raw/*.goalc.txt FILES TO BEASTIE FOR UPDATING');

%% 3) Generate new *goalc.txt files on beastie02 using ktrecon_write_goalc.m
%    (requires reconFrame) and then copy back to local computer


%% Generate gradient_moment_vals.txt / gradient_moment_dirs.txt
%- Create a MATLAB script for each fcmrXXX

%- Manually measure M and S gradient moments
% - Do this manually in the simulator GVE
warning('HAVE YOU MEASURED THE GRADIENT MOMENTS USING GVE?');

%- This is useful as it gives me a record of Vmps and what I did to
% generate the world coordinate version of gradient_moments.
%- Send to beastie02 once done.
warning('YOU MIGHT NEED TO MAKE THE WORLD COORDINATE GRADIENT MOMENT .txt FILES');

