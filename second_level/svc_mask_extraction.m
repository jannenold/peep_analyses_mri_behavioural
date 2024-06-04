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

%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/Frontal_Med_Orb_L.nii,1'
%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/Frontal_Med_Orb_R.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/roi_rvm_tinnermann_0_-32_-44_3mm.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/frontal_midline_bin_de_la_vega.nii,1'
%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/frontal_midline_bin.nii,1'
%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/prgACC.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/Insula_L.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/Insula_R.nii,1'
%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/antIns_L.nii,1'
%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/antIns_R.nii,1'
'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/PAG.nii,1'
%'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/roi_pag_eippert_0_-35_-17_5mm.nii,1'
};
matlabbatch{e}.spm.util.imcalc.output = 'svc_mask';
matlabbatch{e}.spm.util.imcalc.outdir = {'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/'};
matlabbatch{e}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5)>0';
matlabbatch{e}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{e}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{e}.spm.util.imcalc.options.mask = 0;
matlabbatch{e}.spm.util.imcalc.options.interp = 1;
matlabbatch{e}.spm.util.imcalc.options.dtype = 4;

%% Restrict to secondlevel mask space
% e = e +1;
% matlabbatch{e}.spm.util.imcalc.input = {
%                                         '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/svc_mask.nii,1'
%                                         '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/mask_secondlevel.nii,1'
%                                         };
% matlabbatch{e}.spm.util.imcalc.output = 'svc_mask_fov';
% matlabbatch{e}.spm.util.imcalc.outdir = {'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc'};
% matlabbatch{e}.spm.util.imcalc.expression = 'i1.*i2';
% matlabbatch{e}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{e}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{e}.spm.util.imcalc.options.mask = 0;
% matlabbatch{e}.spm.util.imcalc.options.interp = 1;
% matlabbatch{e}.spm.util.imcalc.options.dtype = 4;

%% Smooth Mask
e = e+1;
matlabbatch{e}.spm.spatial.smooth.data = {'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/svc_mask.nii,1'};
matlabbatch{e}.spm.spatial.smooth.fwhm = [1 1 1];
matlabbatch{e}.spm.spatial.smooth.dtype = 0;
matlabbatch{e}.spm.spatial.smooth.im = 0;
matlabbatch{e}.spm.spatial.smooth.prefix = 's_';

%% Create Atlas (with smoothed individual regions so it resembles overal smoothing)
e = e+1;
matlabbatch{e}.spm.util.imcalc.input = {
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/Frontal_Med_Orb_L.nii,1'
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/Frontal_Med_Orb_R.nii,1'
    '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/roi_rvm_tinnermann_0_-32_-44_3mm.nii,1'
    '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/frontal_midline_bin_de_la_vega.nii,1'
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/frontal_midline_bin.nii,1'
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/prgACC.nii,1'
    '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/Insula_L.nii,1'
    '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/Insula_R.nii,1'
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/antIns_L.nii,1'
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/antIns_R.nii,1'
    '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/PAG.nii,1'
    %'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/roi_pag_eippert_0_-35_-17_5mm.nii,1'
    };
matlabbatch{e}.spm.util.imcalc.output = 'svc_mask_atlas';
matlabbatch{e}.spm.util.imcalc.outdir = {'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/heat_pressure/'};
matlabbatch{e}.spm.util.imcalc.expression = '((i1>0).*1+(i2>0).*2+(i3>0).*3+(i4>0).*4+(i5>0).*5)';
matlabbatch{e}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{e}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{e}.spm.util.imcalc.options.mask = 0;
matlabbatch{e}.spm.util.imcalc.options.interp = 1;
matlabbatch{e}.spm.util.imcalc.options.dtype = 4;


spm_jobman('run', matlabbatch);
cd('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/svc/')