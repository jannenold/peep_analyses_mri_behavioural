function [nuissance, n_nuis] = get_move_regressors(func_dir, remain_volumes)

% load the preprocessed movement regressors and return the way used in 
% SPM. One subject at the time.

% The movement slices are already cropped in the same length as my
% runs/onsets. Definition of remaining volumes is done in
% get_nr_of_volumes. remain volumes input = without the unexplained resting
% variance scans at the end of the run!

% Adapted from the PhysIO toolbox (Kasper) included in TAPAS toolbox (Frässle)

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
fm_files    = spm_select('FPList', func_dir, '^rp.*\.txt'); 
nSess       = 2;
n_runs      = 4;

%% motion 24 model

% loop over runs, get the files, compute derivatives and squares
% concat for run length and 

for s = 1:nSess

    for r = 1:n_runs
 
    fm_file      = load(spm_select('FPList', func_dir, sprintf('^rp_asub.*_ses-%02d_task-peep_run-%02d_bold.txt$',s,r)));
    
    % crop the file
    fm_file      = fm_file(1:remain_volumes(r),:);


    % norm the file (for better comparison with other betas)
    SD           = std(fm_file,1);
    mu           = mean(fm_file, 1);
    normed_vals  = (fm_file-mu) ./ SD;

    % square the normed values
    squared_vals = normed_vals.^2;

    % take derivative of values, zeros as first value
    % difference between (t) - (t-1) = slope
    derivs = diff(normed_vals);
    derivs = [zeros(1,size(derivs,2)); derivs];

    % mean correct the derivatives 
    SD           = std(derivs);
    mu           = mean(derivs, 1);
    cor_derivs   = (derivs-mu)./ SD;

    % square the mean corrected derivs
    %squared_derivs = cor_derivs .^2;
    %movements{r} =   [normed_vals, squared_vals, derivs, squared_derivs];

    % only 6 regressors
    movements{r} =   fm_file;

    end
end

nuissance = cat(1, movements{:});
n_nuis    = size(nuissance, 2);

end