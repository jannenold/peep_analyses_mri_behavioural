function complete_FFX_level_peep

% ---------------------------------
% Running complete FFX for subjects
% Model specification:
% Model 1 = HRF
% Model 2 = FIR
% ---------------------------------

clear all; clc;
addpath('/common/apps/spm12-7771');
addpath(genpath('/home/nold/Desktop/PEEP/fMRI_new/'));

%% Specific Settings
subIDs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48];  
exclude_subs    = [8];
session         = [1,2];
model           = 1; % 1 = HRF, 2 = FIR
skernel         = 6;% smoothing kernel for con-images 4mm

%% Select modelname
%modelname = 'complete_onset_stimulus'; %full length, onset start stimulus, end plateau end
%modelname = 'complete_onset_stimulus_condition_contrast'; % 24 condition contrasts (see contrast list)
%modelname = 'complete_onset_plateau';  %full length, onset start plateau
%modelname = 'complete_different_onsets_heat_pressure'; % full length
modelname = 'split_onset_stimulus'; % early and late pain, onset stim start


%% Run FFX Model 
for sub = 1:length(subIDs)
    FFX_peep(subIDs(sub),session,model,modelname);
end

%% VASA
% for sub = 1:length(subIDs)
%     cb_vasa(subIDs(sub),model,modelname);
% end

%% Normalise
for sub = 1:length(subIDs)
    normalise(subIDs(sub),model,modelname);
end

%% Smooth
for sub = 1:length(subIDs)
    smooth(subIDs(sub),model,skernel,modelname);
end
 
%% Warp Masks and prepare for 2nd level
warp_masks(subIDs,model,modelname);
compute_2ndlevel_mask(model,subIDs,exclude_subs,modelname);
%check_masks(subIDs,model,modelname);

end