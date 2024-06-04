function complete_RFX_level_peep
%Compute secondlevel analysis including contrasts for both GLMs
%This includes the following steps:
%
% 1) specify model, i.e. write SPM.mat
% 2) run the model
% 3) specify contrasts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO DO:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; %close all;
addpath('/common/apps/spm12-7771');
addpath(genpath('/home/nold/Desktop/PEEP/fMRI_new/'));

subIDs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 
exclude_subs    = [8];
model           = 1; % 1 = HRF, 2 = FIR;
skern           = 6; % smoothing kernel from FFX

%% Select modelname
%modelname       = ['complete_onset_stimulus_' num2str(skern) 'mm']; %full length, onset start stimulus
%modelname      = 'complete_onset_stimulus_condition_contrast'; % 24 condition contrasts (see contrast list)
%modelname      = 'complete_onset_plateau';  %full length, onset start plateau
%modelname      = 'complete_different_onsets_heat_pressure'; % full length
%modelname      = 'split_onset_stimulus'; % early and late pain, onset stim start
modelname      = 'split_onset_stimulus_fact'; % early and late pain, onset stim start (Factorial)

write_RFX_model(subIDs,exclude_subs,model,modelname,skern);
run_RFX_model(model,modelname,skern);
%compute_contrasts(subIDs,exclude_subs,model,modelname);

end


