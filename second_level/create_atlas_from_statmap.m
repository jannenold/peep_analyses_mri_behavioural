%% ---------------------------------------------
%% Creating SVC Mask
%%----------------------------------------------
clear all
% Input Structures:
%
% PAG: Antaomical from brainstem nav
% antIns_L/antIns_R: left/right anterior insula from Faillenot et al. 2017 comprised of:
%           - insula anterior short gyrus
%           - insula anterior inferior cortex
%           - inusla anterior short gyrus
% RVM: ROI 3mm from Tinnermann et al. 20xx
% Frontal Midline comprised of:
%           - vmPFC (AAL3)
%           - rACC ROI 15 mm from Eippert et al 2009
%
% Thresholded to > 0 (binarised
% Restricted to secondlevel mask
% Smoothing with 2mm Kernel
% -----------------------------------------------------------

%% (Optional) Extract Regions from mask with certain intensitiy labels
% % Load the image
% image = spm_vol('/home/nold/Desktop/PEEP/fMRI_new/masks/neurosynth-mfc-master/neurosynth-mfc-master/images/cluster_labels_k12.nii');
% image_data = spm_read_vols(image);
%
% % Define the intensity label you want to extract
% intensity_label = 4; % Example intensity label
% tolerance = 0.1; % Tolerance range for intensity matching
%
% % Create a binary mask
% binary_mask = abs(image_data - intensity_label) <= tolerance;
%
% % Save the binary mask (optional)
% binary_mask_image = image;
% binary_mask_image.fname = 'vmPFC.nii';
% spm_write_vol(binary_mask_image, binary_mask);

%% Create Mask
e = 1;
matlabbatch{e}.spm.util.imcalc.input = {
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/inter_int_mod_uncorr.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/conj_heat_pressure_int_uncorr.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/neg_inter_int_mod_uncorr.nii,1'
% '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/heat_pressure_sal_fwe.nii,1'
% '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/conj_heat_pressure_fwe.nii,1'
% '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/pressure_heat_sal_fwe.nii,1'

};
matlabbatch{e}.spm.util.imcalc.output = 'comb_heat_pressure_int';
matlabbatch{e}.spm.util.imcalc.outdir = {'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/'};
matlabbatch{e}.spm.util.imcalc.expression = '((i1>0).*1+(i2>0).*2+(i3>0).*4)';
matlabbatch{e}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{e}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{e}.spm.util.imcalc.options.mask = 0;
matlabbatch{e}.spm.util.imcalc.options.interp = 1;
matlabbatch{e}.spm.util.imcalc.options.dtype = 4;

cd('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset/');
spm_jobman('run', matlabbatch);
