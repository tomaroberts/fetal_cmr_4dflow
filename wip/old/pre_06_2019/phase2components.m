function [Pcomp, compFractions_world, compFractions_xyz] = phase2components(P,Vmps,goalcPath)

%% phase2components_world   GENERATE PHASE COMPONENT IMAGES IN WORLD (RAI) COORDINATES
%
% INPUT:
%       - P                     Phase image
%       - Vmps                  First moments in MPS coordinates
%       - goalcPath             Path to goalc.txt file
%
% OUTPUT:
%       - Pcomp                 Component phase images
%       - compFractions_world   Vworld as fractions
%       - compFractions_xyz     Vxyz as fractions
%




%% getgoalc parameters
gc = get_pcmr_orientation_parameters( goalcPath );


%% First moment values from GVE
% Vm = 9.95;
% Vp = 0;
% Vs = -10.98;
% 
% Vmps = [Vm; Vp; Vs];


%% Convert First Moments to xyz/World Coordinates
[Vworld, Vxyz] = vmps2vworld(Vmps,gc);


%% Component phase images

momNorm = @(x) x ./ norm(x); %calculate components of vector norm of gradient first moments

compFractions_world = momNorm(Vworld);
for ii = 1:3
    Pcomp(:,:,:,ii,1) = P .* compFractions_world(ii)';
end

compFractions_xyz = momNorm(Vxyz);
for ii = 1:3
    Pcomp(:,:,:,ii,2) = P .* compFractions_xyz(ii)';
end


% momfrac = @(x) abs(x ./ sum(abs(x))); %calculate moments as percentages

% compFractions_world = momfrac(Vworld);
% for ii = 1:3
%     Pcomp(:,:,:,ii,1) = P .* compFractions_world(ii)';
% end
% 
% compFractions_xyz = momfrac(Vxyz);
% for ii = 1:3
%     Pcomp(:,:,:,ii,2) = P .* compFractions_xyz(ii)';
% end

% %% view
% s = 43;
% cscale = [-pi/8,pi/8];
% 
% imtar(P(:,:,s),cscale(1),cscale(2)); title('original');
% 
% im1 = squeeze([Pcomp(:,:,s,1), Pcomp(:,:,s,2), Pcomp(:,:,s,3)]);
% imtar(im1,cscale(1),cscale(2));
% title('components (RL-AP-FH)');




end %end fn