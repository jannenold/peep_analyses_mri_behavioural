%for PAG: 
% '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/rois/brainstem_nav/brainstem_nav.nii'
%for ACC, Insula, Operc, vmPFC; 
% '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/rois/aal3_atlas/AAL3/AAL3v1_1mm.nii'

xA=spm_atlas('load','/home/nold/Desktop/PEEP/fMRI_new/masks/aal/AAL3v1_for_SPM12/AAL3/AAL3v1_1mm.nii');
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN_old/derivatives_JN/spm_secondlevel/rois/brainstem_nav/brainstem_nav.nii'); 

S=spm_atlas('select',xA);

for i = 1:size(S,2)
fname=strcat(S{i},'.nii');
VM=spm_atlas('mask',xA,S{i});
VM.fname=fname;
spm_write_vol(VM,spm_read_vols(VM));
end

