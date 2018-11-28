function [Vworld, Vxyz] = vmps2vworld(Vmps,gc)

%% vmps2vworld   CONVERT GRADIENT MOMENT VECTORS FROM MPS TO WORLD COORDINATES
%
% INPUT:
%       - Vmps          First moments measured in GVE
%       - gc            goalc structure containing m_/p_/s_/orient
%
% OUTPUT:
%       - Vworld        First moments in world coordinates
%       - Vxyz          First moments in scanner (xyz) coordinates
%
% TO DO:
% - eventually change so that it can calculate Vmps from goalc structure
% - need to establish whether gradient dirs are consistent with m_orient/s_orient
% - is Cprime correct? Should it be RAS (ie: transformed to same as .nii?)


%% First moments
Vm = Vmps(1);
Vp = Vmps(2);
Vs = Vmps(3);

Vmps = [Vm; Vp; Vs];


%% F --- First moment orientations
% think I can use this to automatically determine gradient moment directions
F.m = [gc.mc0_str, gc.m0_str ];  % should be [-ve, +ve]
F.s = [gc.s_ex_str, gc.r0_str ]; % should be [+ve, -ve]

if sign(F.m) ~= [-1, 1]
    warning('Readout gradients are orientated in an expected way!');
end

if sign(F.s) ~= [1, -1]
    warning('Slice-select gradients are orientated in an unexpected way!');

    % update Vmps if slice-select is flipped
    if sign(F.s) == [-1, 1]
        Vs = Vs * -1;
        Vmps = [Vm; Vp; Vs];
        warning('Slice-select gradient flipped: Vs = -Vs');
    end
end


%% O_MATRIX --- converts direction of slice from MPS to xyz coordinates
%- get from goalc: 'O_MATRIX `locations[0]'
%- eg: for transverse slice:
%- O_MATRIX = [m_orient; p_orient; s_orient];
%- m_orient = [0 -1 0]; <- means M in -y direction (= LR in Philips)
%- p_orient = [1 0 0];  <- means P in +x direction (= PA in Philips)
%- s_orient = [0 0 -1]; <- means S in -z direction (= HF in Philips)
m_orient = gc.m_orient;
p_orient = gc.p_orient;
s_orient = gc.s_orient;
O = [m_orient; p_orient; s_orient];


%% Define Gradient Moment direction
%- relationship between direction of gradients and +ve/-ve velocity encoding
%- Either positive or negative
%- Worked out that needs to be negative by matching to equivalent QFlow
D = [-1, 0, 0; 0 -1, 0; 0 0 -1];


%% First moments in xyz coordinates
Vxyz = D * O' * Vmps;


%% Transform from scanner to world
%- Philips defines scanner as:
%- +x = PA
%- +y = RL
%- +z = FH (IS)

%- World/Patient definition (ie: RAI):
%%%%% I think this might need to match .nii coord system? ie: RAS?
%- +x` = RL
%- +y` = AP
%- +z` = FH (IS)

%- conversion from Philips to normal world:
%- scanner : world
%-       x : -y`
%-       y :  x`
%-       z :  z`

Cprime = [0 -1 0; 1 0 0; 0 0 1];


%% Perform transform to world coordinates
%- currently, based on Cprime above:
%- Vworld(1) = RL
%- Vworld(2) = AP
%- Vworld(3) = FH
Vworld = Cprime * Vxyz;


% end vmps2vworld(...).m
end