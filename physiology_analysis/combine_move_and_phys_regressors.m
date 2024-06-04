function combine_move_and_phys_regressors(sub,exclude,nSess,n_runs,do_bs)

physio_dir = '/projects/crunchie/nold/PEEP/physiology/MAIN/noise_regressors';

% exclude relevant subjects
if ~isempty(exclude)
    sub = sub(~ismember(sub,exclude));
end

for g=1:numel(sub)
    l = 1;
    clear n_nuis

    for s = 1:nSess
        for r = 1:n_runs

            clear nui_reg motion physio

            if do_bs
                filename     = sprintf('BS_complete_reg_sub-%02d_ses_%02d_run_%02.2d.txt',sub(g),s,r);
                filename_mat     = sprintf('BS_complete_reg_sub-%02d_ses_%02d_run_%02.2d.mat',sub(g),s,r);
            else
                filename     = sprintf('complete_reg_sub-%02d_ses_%02d_run_%02.2d.txt',sub(g),s,r);
                filename_mat     = sprintf('complete_reg_sub-%02d_ses_%02d_run_%02.2d.mat',sub(g),s,r);
            end

            % ---------------------------------------------------------------------
            %               Motion
            % ---------------------------------------------------------------------

            if do_bs
                motion_file_name = sprintf('sub-%02d_BS_motion_ses-%02d_run-%02d.mat',sub(g),s,r);
            else
                motion_file_name = sprintf('sub-%02d_motion_ses-%02d_run-%02d.mat',sub(g),s,r);
            end


            if exist([physio_dir filesep motion_file_name])
                load([physio_dir filesep motion_file_name]);
                dim_motion = size(motion,1); % check number of scans
            else
                fprintf('No Motion File for sub-%02d session %02d run %02d',sub(g),s,r);
            end

            % ---------------------------------------------------------------------
            %               spikes
            % ---------------------------------------------------------------------

            spike_file_name = sprintf('sub-%02d_nui_reg_mov_ses-%02d_run-%02d.mat',sub(g),s,r);

            if exist([physio_dir filesep spike_file_name])
                load([physio_dir filesep spike_file_name]);
                dim_nui_reg = size(nui_reg,1); % check number of scans
            else
                fprintf('\nNo Spikes File for sub-%02d session-%02d run-%02d\n',sub(g),s,r);
                nui_reg = [];
            end

            % ---------------------------------------------------------------------
            %                 Physio regressors
            % ---------------------------------------------------------------------

            physio_file_name = sprintf('sub-%02d_physio_ses-%02d_run-%02d.mat',sub(g),s,r);

            if exist([physio_dir filesep physio_file_name])
                load([physio_dir filesep physio_file_name]);
                dim_physio = size(physio,1); % check number of scans
            else
                fprintf('\nNo Physio File for sub-%02d session-%02d run-%02d\n',sub(g),s,r);
                physio = [];
            end

            % ---------------------------------------------------------------------
            %                 add all nuissance regressors
            % ---------------------------------------------------------------------
            complete_reg =   [motion,physio,nui_reg];

            % save as .txt
            save(fullfile(physio_dir,filename),'complete_reg','-ascii');

            % also save as .mat 
            save(fullfile(physio_dir,filename_mat),'complete_reg');
            
            % save number of nuissance regressors as mat file for later
            % contrast estimation  
            n_nuis(l) = size(complete_reg,2);
            l = l+ 1;
        end
    end
    save(fullfile(physio_dir,sprintf('n_nuis_sub-%02d',sub(g))),'n_nuis');

end

