function BS_complete_RFX_level_peep
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
addpath(genpath('/home/nold/Desktop/PEEP/fMRI/'));

subIDs          = [1:4 6:9 11 13:15 17:19 21:32 34:36 38:41 43:48]; 
exclude_subs     = [4,8,18,22,29,32,41,45];
model           = 1; % 1 = HRF, 2 = FIR;
skern           = 4; % smoothing kernel from FFX

%% Select modelname
modelname       = ['complete_onset_stimulus_' num2str(skern) 'mm']; %full length, onset start stimulus
%modelname      = 'complete_onset_stimulus_condition_contrast'; % 24 condition contrasts (see contrast list)
%modelname      = 'complete_onset_plateau';  %full length, onset start plateau
%modelname      = 'complete_different_onsets_heat_pressure'; % full length
%modelname      = 'split_onset_stimulus'; % early and late pain, onset stim start


%BS_write_RFX_model(subIDs,exclude_subs,model,modelname,skern);
BS_run_RFX_model(model,modelname,skern);


end


