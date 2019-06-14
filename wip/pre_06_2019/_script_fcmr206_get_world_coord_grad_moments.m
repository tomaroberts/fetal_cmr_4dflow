% fcmr206 --- creating gradient_moment .txt files
%- not entirely sure that I have them in correct world space orientation
%  relative to the .nii of the stack

cd('C:\Users\tr17\Documents\Projects\PC_Fetal_CMR\Data\fcmr206\ktrecon');

gcFiles = {'s17_goalc','s18_goalc','s20_goalc','s21_goalc'};

for ii = 1:numel(gcFiles)
    gc = get_pcmr_orientation_parameters( [gcFiles{ii} '_TAR.txt'] ); %_TAR to differentiate from Josh's original goalc files.
    
    m_orient = gc.m_orient;
    p_orient = gc.p_orient;
    s_orient = gc.s_orient;
    O(:,:,ii) = [m_orient; p_orient; s_orient];
    
    [Vworld(:,ii), Vxyz(:,ii)] = vmps2vworld(Vmps,gc);
end


% grad_moment directions
for ii = 1:size(Vworld,2)
    Vworld_unit(:,ii) = Vworld(:,ii)./norm(Vworld(:,ii));
end

% grad_moment magnitude
for ii = 1:size(Vworld,2)
    Vworld_value(:,ii) = norm(Vworld(:,ii));
end

% save to .txt files compatible with SVRTK
cd ../vel_vol

fileID = fopen('grad_moment_dirs.txt','w');
fprintf(fileID,'%.4f\t',Vworld_unit(:));
fclose(fileID);

fileID = fopen('grad_moment_vals.txt','w');
fprintf(fileID,'%.4f\t',Vworld_value(:));
fclose(fileID);