function BS_complete_FFX_level_peep

% ---------------------------------
% Running complete FFX for subjects
% Model specification:
% Model 1 = HRF
% Model 2 = FIR
% ---------------------------------
%% TO DO:
% - include VASA after FFX (CB)
%-----------------------------------

clear all; clc;
addpath('/common/apps/spm12-7771');
addpath(genpath('/home/nold/Desktop/PEEP/fMRI/'));

%% Specific Settings
subIDs          = [3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48];  
exclude_subs     = [4,8,18,22,29,32,41,45];
session         = [1,2];
model           = 1; % 1 = HRF, 2 = FIR
skernel         = 4;% smoothing kernel for con-images 4mm

%% Select modelname
modelname = 'complete_onset_stimulus'; %full length, onset start stimulus, end plateau end
%modelname = 'complete_onset_stimulus_condition_contrast'; % 24 condition contrasts (see contrast list)
%modelname = 'complete_onset_plateau';  %full length, onset start plateau
%modelname = 'complete_different_onsets_heat_pressure'; % full length
%modelname = 'split_onset_stimulus'; % early and late pain, onset stim start


%% Run FFX Model 
%for sub = 1:length(subIDs)
%    BS_FFX_peep(subIDs(sub),exclude_subs,session,model,modelname);
%end


%% Normalise
%BS_normalise(subIDs,exclude_subs,model,modelname);


%% Smooth
BS_smooth(subIDs,exclude_subs,model,skernel,modelname);


end