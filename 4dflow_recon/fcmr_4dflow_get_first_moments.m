function fcmr_4dflow_get_first_moments( fcmrDir, varargin )
%FCMR_4DFLOW_GET_FIRST_MOMENTS  get and save gradient first moments
%
%   FCMR_4DFLOW_GET_FIRST_MOMENTS( fcmrDir, fcmrNum, 'param', val )
%       Get gradient first moments from .raw files and save to text files
%       for use in SVRTK reconstruction
%       - IMPORTANT: currently requires manual measurement of gradient
%       first moments using GVE in Philips simulator.
% 
%   Requires:
%       - .raw data in ...fcmrDir/raw
%       - Measurement of gradient first moments (i.e: using GVE)
%       
%   Input:
%       fcmrDir             path to fcmrXXX directory
%
%   Optional Parameter-value Pairs:
%       Vmps                gradient first moments in mps (Measurement/Phase/Slice) coordinates (Philips terminology)
%
%   Output:
%       Various .text files containing reorientated gradient first moments
%
%   Example usage:
%       - Measure gradient first moments in GVE
%       fcmr_4dflow_postprocessing( 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr194' )
%
%   TODO:
%       - automatically calculate gradient first moments from .raw files,
%       rather than require manual measurement in GVE
%
%   See also:
%       

% Tom Roberts (t.roberts@kcl.ac.uk)


%% Parse Input

default.Vmps              = [5.96, 0, -4.35]; % default values for iFIND protocol

p = inputParser;

if  verLessThan('matlab','8.2')
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParamValue( parseobj, argname, defaultval, validator );
else
    add_param_fn = @( parseobj, argname, defaultval, validator ) addParameter( parseobj, argname, defaultval, validator );
end

addRequired(  p, 'fcmrDir' );

add_param_fn( p, 'Vmps', default.Vmps, ...
        @(x) validateattributes( x, {'double'}, ...
        {}, mfilename ) );

parse( p, fcmrDir, varargin{:} );

Vmps               = p.Results.Vmps;


%% Directory admin
% fcmrDir = 'C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr201';
rawDir = '/raw';
ktreconDir = '/ktrecon';
dataDir = '/data';


%% First moment values from GVE
%- IMPORTANT: MEASURE THESE MANUALLY USING GVE IN SIMULATOR
Vm = Vmps(1);
Vp = Vmps(2);
Vs = Vmps(3);

Vmps = [Vm; Vp; Vs];

if isempty(Vm) || isempty(Vs)
    error('Vm / Vs empty. You must supply gradient moment values measured in GVE.');
end


%% Rename goalc.txt files
cd([fcmrDir rawDir]);
gcNames = dir('*goalc.txt');
for ss = 1:numel(gcNames)
    sIDs{ss} = gcNames(ss).name(21:22);
    copyfile( gcNames(ss).name , [fcmrDir ktreconDir '/s' sIDs{ss} '_goalc_TAR.txt' ] ); %_TAR to differentiate from Josh's original goalc files
end
clear gcNames


%% Convert to world/xyz coordinates
cd([fcmrDir ktreconDir]);

gcFiles = dir('*goalc_TAR.txt');

for ii = 1:numel(gcFiles)
    gc = get_pcmr_orientation_parameters( gcFiles(ii).name );
    
    % get orient objects for record
    m_orient = gc.m_orient;
    p_orient = gc.p_orient;
    s_orient = gc.s_orient;
    O(:,:,ii) = [m_orient; p_orient; s_orient];
    
    % get M / S gradient strengths for record
    gr_str(ii).mc0_str  = gc.mc0_str;
    gr_str(ii).m0_str   = gc.m0_str;
    gr_str(ii).m3_str   = gc.m3_str;
    gr_str(ii).s_ex_str = gc.s_ex_str;
    gr_str(ii).r0_str   = gc.r0_str;
    
    [Vworld(:,ii), Vxyz(:,ii)] = vmps2vworld(Vmps,gc);
end


% grad_moment directions as unit vector
for ii = 1:size(Vworld,2)
    Vworld_unit(:,ii) = Vworld(:,ii)./norm(Vworld(:,ii));
end

% grad_moment magnitude
for ii = 1:size(Vworld,2)
    Vworld_value(:,ii) = norm(Vworld(:,ii));
end


%% save grad_moment.txt files compatible with SVRTK
%  (and ORIENT object)
cd ../data

% For SVRTK
fileID = fopen('grad_moment_dirs.txt','w');
fprintf(fileID,'%.4f\t',Vworld_unit(:));
fclose(fileID);

fileID = fopen('grad_moment_vals.txt','w');
fprintf(fileID,'%.4f\t',Vworld_value(:));
fclose(fileID);

% For reference
fileID = fopen('grad_moments_measurementCoords.txt','w');
for ii = 1:numel(gcFiles)
    fprintf(fileID,['== s' sIDs{ii} ' ==\n']);
    fprintf(fileID,'%.4f %.4f %.4f\n', Vmps');
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('grad_moments_worldCoords.txt','w');
for ii = 1:numel(gcFiles)
    fprintf(fileID,['== s' sIDs{ii} ' ==\n']);
    fprintf(fileID,'%.4f %.4f %.4f\n', Vworld(:,ii)');
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('ORIENT.txt','w');
for ii = 1:numel(gcFiles)
    fprintf(fileID,['== s' sIDs{ii} ' ==\n']);
    fprintf(fileID,'%.4f %.4f %.4f\n', O(:,:,ii));
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('gr_str.txt','w');
for ii = 1:numel(gcFiles)
    fprintf(fileID,['== s' sIDs{ii} ' ==\n']);
    fprintf(fileID,['mc0_str = %.4f\n'], gr_str(ii).mc0_str);
    fprintf(fileID,['m0_str = %.4f\n'], gr_str(ii).m0_str);
    fprintf(fileID,['m3_str = %.4f\n'], gr_str(ii).m3_str);
    fprintf(fileID,['s_ex_str = %.4f\n'], gr_str(ii).s_ex_str);
    fprintf(fileID,['r0_str = %.4f\n'], gr_str(ii).r0_str);
    fprintf(fileID,'\n');
end
fclose(fileID);


end %fcmr_4dflow_get_first_moments(...)



