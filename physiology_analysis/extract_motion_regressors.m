function extract_motion_regressors(subs,exclude,n_runs,nSess,do_bs)

fprintf('\n\n--------------------------------------------------------\n');
fprintf('                       Motion                               \n');
fprintf('------------------------------------------------------------\n');
outdir = '/projects/crunchie/nold/PEEP/physiology/MAIN/noise_regressors';

if ~isempty(exclude)
    subs = subs(~ismember(subs,exclude));
end

for g=1:numel(subs)
    for s = 1:nSess
        for r = 1:n_runs

            sub_dir = sprintf('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_preprocessing/sub-%02d/func',subs(g));


            fprintf('Extracting Motion Parameters sub-%02d ses-%02d run%02d\n',subs(g),s,r);
            fprintf('=================================================\n');

            % ---------------------------------------------------------------------
            %               movement
            % ---------------------------------------------------------------------
            if do_bs
                fm_file      = load(spm_select('FPList', sub_dir, sprintf('^rp_basub.*_ses-%02d_task-peep_run-%02d_bold.txt$',s,r)));
            else
                fm_file      = load(spm_select('FPList', sub_dir, sprintf('^rp_asub.*_ses-%02d_task-peep_run-%02d_bold.txt$',s,r)));
            end
            % crop the file
            %fm_file      = fm_file(1:remain_volumes(r),:);

            % norm the file (for better comparison with other betas)
            SD           = std(fm_file,1);
            mu           = mean(fm_file, 1);
            normed_vals  = (fm_file-mu) ./ SD;

            % square the normed values
            squared_vals = normed_vals.^2;

            % take derivative of values, zeros as first value
            %difference between (t) - (t-1) = slope
            derivs = diff(normed_vals);
            derivs = [zeros(1,size(derivs,2)); derivs];

            % mean correct the derivatives
            SD           = std(derivs);
            mu           = mean(derivs, 1);
            cor_derivs   = (derivs-mu)./ SD;

            % square the mean corrected derivs
            squared_derivs = cor_derivs .^2;

            motion =   [normed_vals, squared_vals, derivs, squared_derivs];
            if do_bs
                outFilename = sprintf('sub-%02.0f_BS_motion_ses-%02.0f_run-%02.0f.mat',subs(g),s,r);
            else
                outFilename = sprintf('sub-%02.0f_motion_ses-%02.0f_run-%02.0f.mat',subs(g),s,r);
            end

            save(fullfile(outdir,outFilename),'motion');

        end
    end
end