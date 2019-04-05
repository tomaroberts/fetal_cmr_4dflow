%% _investigating_flowphantom6_v_fetal_gradient_moments
%- trying to see that they make sense


%% Flow Phantom 6
sp = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\2018_11_01_Flow_Phantom6\PC-bSSFP\raw_data';
cd(sp);


% Get goalc.txt filepaths:
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

% 5 stacks - tra/sag/cor/tra_obl/sag_obl
selectedStacks = [1,2,3,7,8];
% selectedStacks = 1:9;

% First moment values from GVE
Vm = 9.95;
Vp = 0;
Vs = -10.98;

Vmps_fp6 = [Vm; Vp; Vs];

% Investigate goalc parameters
gc_fp6 = {};
for ii = 1:numel(selectedStacks)
	gc_fp6{ii} = get_pcmr_orientation_parameters( [scanName{selectedStacks(ii)} '_goalc.txt'] );
    
    O_fp6(1,:,ii) = gc_fp6{ii}.m_orient;
    O_fp6(2,:,ii) = gc_fp6{ii}.p_orient;
    O_fp6(3,:,ii) = gc_fp6{ii}.s_orient;
    
    % slice select gradient strength - gives directional information
    s_ex_str_fp6(ii) =  gc_fp6{ii}.s_ex_str;
    r0_str_fp6(ii) =  gc_fp6{ii}.r0_str;
end

xEncStr_fp6 = sum(abs(nonzeros(O_fp6([1,3],1,:))),1);
yEncStr_fp6 = sum(abs(nonzeros(O_fp6([1,3],2,:))),1);
zEncStr_fp6 = sum(abs(nonzeros(O_fp6([1,3],3,:))),1);
EncStr_fp6 = [xEncStr_fp6, yEncStr_fp6, zEncStr_fp6];

% Convert Vmps to world coordinates for each stack
for ii = 1:numel(selectedStacks)
	Vworld_fp6(:,ii) = vmps2vworld(Vmps_fp6, get_pcmr_orientation_parameters( [scanName{selectedStacks(ii)} '_goalc.txt'] ) );
end

% grad_moment directions
for ii = 1:size(Vworld_fp6,2)
    Vworld_fp6_unit(:,ii) = Vworld_fp6(:,ii)./norm(Vworld_fp6(:,ii));
end

% grad_moment magnitude
for ii = 1:size(Vworld_fp6,2)
    Vworld_fp6_value(:,ii) = norm(Vworld_fp6(:,ii));
end

%% fcmr191
%-s16 = tra
%-s17 = sag
%-s18 = cor
%-s19 = 2ch (2-chamber?)
%-s20 = sa (short axis?)

cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr191\ktrecon');

% goalc.txt filenames
gcFiles = {'s16_goalc','s17_goalc','s18_goalc','s19_goalc','s20_goalc'};

% First moment values from GVE
Vm = 5.96;
Vp = 0;
Vs = -4.35;

Vmps = [Vm; Vp; Vs];

% Investigate goalc parameters
gc_fcmr191 = {};
for ii = 1:numel(gcFiles)
	gc_fcmr191{ii} = get_pcmr_orientation_parameters( [gcFiles{ii} '_TAR.txt'] );
    
    O_fcmr191(1,:,ii) = gc_fcmr191{ii}.m_orient;
    O_fcmr191(2,:,ii) = gc_fcmr191{ii}.p_orient;
    O_fcmr191(3,:,ii) = gc_fcmr191{ii}.s_orient;   
    
    % slice select gradient strength - gives directional information
    s_ex_str_fcmr191(ii) =  gc_fcmr191{ii}.s_ex_str;
    r0_str_fcmr191(ii) =  gc_fcmr191{ii}.r0_str;
end

xEncStr_fcmr191 = sum(abs(nonzeros(O_fcmr191([1,3],1,:))),1);
yEncStr_fcmr191 = sum(abs(nonzeros(O_fcmr191([1,3],2,:))),1);
zEncStr_fcmr191 = sum(abs(nonzeros(O_fcmr191([1,3],3,:))),1);
EncStr_fcmr191 = [xEncStr_fcmr191, yEncStr_fcmr191, zEncStr_fcmr191]

% Convert Vmps to world coordinates for each stack
for ii = 1:numel(gcFiles)
    Vworld_fcmr191(:,ii) = vmps2vworld(Vmps, get_pcmr_orientation_parameters( [gcFiles{ii} '_TAR.txt'] ) );
end

% grad_moment directions
for ii = 1:size(Vworld_fcmr191,2)
    Vworld_fcmr191_unit(:,ii) = Vworld_fcmr191(:,ii)./norm(Vworld_fcmr191(:,ii));
end

% grad_moment magnitude
for ii = 1:size(Vworld_fcmr191,2)
    Vworld_fcmr191_value(:,ii) = norm(Vworld_fcmr191(:,ii));
end