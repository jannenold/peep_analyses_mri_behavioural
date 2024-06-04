%% Extract physiology and movement (spikes)

% Extracts the following regressors and stores them in 
% /projects/crunchie/nold/PEEP/physiology/MAIN/noise_regressors
% 
% - Physiology: Pulse and Respiration
% - Spike Movement: Regressor file for spike movements (only for subs where
% it applies
% - Motion: 24 motion regressors (6 normed, 6 squared, 6 derivs, 6 squared derivs)
% 

% The 6 realignment parameters can be augmented by their derivatives (in
% total 12 parameters) + the squares of parameters and derivatives (in
% total 24 parameters), leading to a Volterra expansion as in
%
%     Friston KJ, Williams S, Howard R, Frackowiak
%     RS, Turner R. Movement-related effects in fMRI
%     time-series. Magn Reson Med. 1996;35:346?355.)

% derivs: temporal difference of realignment parameters,
%         prepending one line of zeros (for 1st scan)

% The physio regressors should also already have the correct length!

% -------------------
% TO DO:
% -------------------


clear all; clc;
addpath('/common/apps/spm12-7771');
addpath(genpath('/home/nold/Desktop/PEEP/fMRI_new/'));

subIDs = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48];
subIDs = 8
exclude_subs = [];
%exclude_subs = [4,8,18,22,29,32,41,45];
nRun = 4;
nSess = 2;
session = [1,2];

do_spikes       = 1; % extract spike movements
do_physio       = 0; % extract physio (pulse,resp)
do_motion       = 0; % extract 24 motion params
combine_regs    = 0; % combine regressors to txt file for later analyses
do_bs           = 0; % set to zero to choose cortical movement reg

% -----------------------------------
% Spike Movement Regressor
% -----------------------------------
if do_spikes
    inspect_movement = 1; % should always be set to 1
    plot_movement = 0;
    generate_reg = 1;
    threshSpike = 0.6; % mm threshold for spike detection
    threshSess = 3;    % mm threshold fo between run/session movement
    inspect_movement_BIDS(subIDs,nRun,nSess,inspect_movement,plot_movement,generate_reg,threshSpike,threshSess);
end

% -----------------------------------
% Physio Files Regressor (Pulse, Resp)
% -----------------------------------
if do_physio
    nScans = 400;
    nDummy = 4;
    extract_physio(subIDs,exclude_subs,nRun,session,nScans,nDummy);
end

%------------------------------------
% 24 Motion Regressors
%------------------------------------
if do_motion
    extract_motion_regressors(subIDs,exclude_subs,nRun,nSess,do_bs);
end

%-------------------------------------
% Combine all regressors and save as txt
%-------------------------------------
if combine_regs
    combine_move_and_phys_regressors(subIDs,exclude_subs,nSess,nRun,do_bs);
end