
%% studyPath
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


%% getgoalc parameters for each scan

tra.p.gc = get_pcmr_orientation_parameters( [scanName{1} '_goalc.txt'] );
sag.p.gc = get_pcmr_orientation_parameters( [scanName{2} '_goalc.txt'] );
cor.p.gc = get_pcmr_orientation_parameters( [scanName{3} '_goalc.txt'] );

tra.m.gc = get_pcmr_orientation_parameters( [scanName{4} '_goalc.txt'] );
sag.m.gc = get_pcmr_orientation_parameters( [scanName{5} '_goalc.txt'] );
cor.m.gc = get_pcmr_orientation_parameters( [scanName{6} '_goalc.txt'] );

tra.o.gc = get_pcmr_orientation_parameters( [scanName{7} '_goalc.txt'] );
sag.o.gc = get_pcmr_orientation_parameters( [scanName{8} '_goalc.txt'] );
cor.o.gc = get_pcmr_orientation_parameters( [scanName{9} '_goalc.txt'] );


%% First moment values
%- currently measured in GVE
%- TO DO: ideally automate this from goalc file by working out first moments of MPS 
Vm = 6.21;
Vp = 0;
Vs = -6.65;

Vmps = [Vm; Vp; Vs];


%% Convert First Moments to World Coordinates
[tra.p.Vworld, tra.p.Vxyz] = vmps2vworld(Vmps,tra.p.gc);
[sag.p.Vworld, sag.p.Vxyz] = vmps2vworld(Vmps,sag.p.gc);
[cor.p.Vworld, cor.p.Vxyz] = vmps2vworld(Vmps,cor.p.gc);

[tra.m.Vworld, tra.m.Vxyz] = vmps2vworld(Vmps,tra.m.gc);
[sag.m.Vworld, sag.m.Vxyz] = vmps2vworld(Vmps,sag.m.gc);
[cor.m.Vworld, cor.m.Vxyz] = vmps2vworld(Vmps,cor.m.gc);

[tra.o.Vworld, tra.o.Vxyz] = vmps2vworld(Vmps,tra.o.gc);
[sag.o.Vworld, sag.o.Vxyz] = vmps2vworld(Vmps,sag.o.gc);
[cor.o.Vworld, cor.o.Vxyz] = vmps2vworld(Vmps,cor.o.gc);







