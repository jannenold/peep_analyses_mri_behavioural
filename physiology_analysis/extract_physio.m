function extract_physio(subs,exclude,n_runs,num_sess,n_scans,dummy)

fprintf('\n\n--------------------------------------------------------\n');
fprintf('                       Physio                               \n');
fprintf('------------------------------------------------------------\n');

if ~isempty(exclude)
    subs = subs(~ismember(subs,exclude));
end


basedir = '/projects/crunchie/nold/PEEP/physiology/MAIN';
preprocdir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_preprocessing';

for g = 1:size(subs,2)
    name = sprintf('sub0%02.2d',subs(g));


    for se = 1:numel(num_sess)

   
        for r = 1:n_runs

            fprintf('\n');
            fprintf('Extracting %s in session ses-%02d for run-%02d\n',name,num_sess(se),r);
            disp('===================================================================');

            physiodir  = [basedir filesep 'original' filesep sprintf('ses-%02d',num_sess(se)) filesep];
            [physio,~,~,~,nscans_ind] = get_physio(physiodir,name,r,n_scans);
            savedir = [basedir filesep 'noise_regressors' filesep];

            % Determine dummy scans in the end (differ between runs since scanner is stopped manually)
            % Load movement file for reference dimensions
            rpFile = sprintf('rp_asub-%02d_ses-%02d_task-peep_run-%02d_bold.txt',subs(g),num_sess(se),r);
            a = load(fullfile(preprocdir,sprintf('sub-%02d',subs(g)), 'func',rpFile));
            dim_mov = size(a,1);

            if isequal(dim_mov+(dummy+1),nscans_ind)
                dummy_end = 1;
            elseif diff([dim_mov+(dummy+1),nscans_ind])==-1
                dummy_end = 0;
            elseif diff([dim_mov+(dummy+1),nscans_ind])==1
                dummy_end = 2;
            elseif diff([dim_mov+(dummy+1),nscans_ind])==2
                dummy_end = 3;
            elseif diff([dim_mov+(dummy+1),nscans_ind])==3
                dummy_end = 4;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==4
                dummy_end = 5;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==5
                dummy_end = 6;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==6
                dummy_end = 7;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==7
                dummy_end = 8;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==8
                dummy_end = 9;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==9
                dummy_end = 10;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==10
                dummy_end = 11;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==11
                dummy_end = 12;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==12
                dummy_end = 13;
                warning('Difference in Physio File Scans!');
            elseif diff([dim_mov+(dummy+1),nscans_ind])==13
                dummy_end = 14;
                warning('Difference in Physio File Scans!');
            else
                error('Check scans in Physio file!');

            end


            if ~isempty(physio)
                physio = physio(dummy+1:nscans_ind-dummy_end,:);
                save(fullfile(savedir,sprintf('sub-%02.0f_physio_ses-%02d_run-%02d.mat',subs(g),num_sess(se),r)),'physio');
            end

        end
    end
end

end

