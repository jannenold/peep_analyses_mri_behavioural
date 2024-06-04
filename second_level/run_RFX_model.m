function run_RFX_model(model,modelname,skern)
%Estimate the pre-written SPM.mat for the secondlevel

%path business
if model == 1
    model_prefix = 'HRF';
elseif model == 2
    model_prefix = 'FIR';
end
estimate_contrasts_heat_pressure= 1;
fact = 1;

% more path business
if strcmp(modelname,'split_onset_stimulus')&& ~estimate_contrasts_heat_pressure
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm']) && ~estimate_contrasts_heat_pressure
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_plateau')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_different_onsets_heat_pressure')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_stimulus_condition_contrast')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm']) && estimate_contrasts_heat_pressure
    modelstr = [model_prefix '_' modelname '_heat_pressure'];
elseif strcmp(modelname,['split_onset_stimulus']) && estimate_contrasts_heat_pressure
    modelstr = [model_prefix '_' modelname '_6mm_heat_pressure'];
elseif strcmp(modelname,['split_onset_stimulus_fact']) && estimate_contrasts_heat_pressure
    modelstr = ['FACT_',model_prefix '_' modelname '_6mm_heat_pressure'];
end



[path]  = get_study_specs;
RFXDir           = path.secondlevelDir;

if model == 1

    if strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm'])
        mats = dir(fullfile(RFXDir,modelstr));
        mats = {mats(3:end).name};
        n_mats    = length(mats);
    elseif strcmp(modelname,'split_onset_stimulus')
        mats = dir(fullfile(RFXDir,modelstr));
        mats = {mats(3:end).name};
        n_mats    = length(mats);
    elseif strcmp(modelname,'complete_onset_stimulus_condition_contrast')
        n_mats    = 1;
    elseif strcmp(modelname,'split_onset_stimulus_fact')
        n_mats    = 1;
    end

    for e = 1:n_mats

        if strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm'])
            outdir = fullfile(RFXDir, modelstr, mats{e});
        elseif strcmp(modelname,'split_onset_stimulus')
            outdir = fullfile(RFXDir, modelstr, mats{e});
        elseif strcmp(modelname,'complete_onset_stimulus_condition_contrast')
            outdir = fullfile(RFXDir, modelstr);
        elseif strcmp(modelname,'split_onset_stimulus_fact')
              outdir = fullfile(RFXDir, [modelstr]);
        end

        spm_mat = fullfile(outdir, 'SPM.mat');
        assert(exist(spm_mat, 'file') > 0, 'No such file: %s', spm_mat);
        matlabbatch{e}.spm.stats.fmri_est.spmmat = cellstr(spm_mat);
        matlabbatch{e}.spm.stats.fmri_est.method.Classical = 1;
        %matlabbatch{e}.spm.stats.con.spmmat = cellstr(spm_mat);
        %matlabbatch{e}.spm.stats.con.consess{1}.tcon.name = 'con';
        %matlabbatch{e}.spm.stats.con.consess{1}.tcon.weights = 1;
        %matlabbatch{e}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        %matlabbatch{e}.spm.stats.con.delete = 0;
   
    end
    

elseif model == 2

    mats = dir(fullfile(RFXDir,modelstr));
    mats = {mats(3:end).name};
    n_mats    = length(mats);
    e = 1;
    outdir = fullfile(RFXDir, modelstr);
    spm_mat = fullfile(outdir, 'SPM.mat');
    assert(exist(spm_mat, 'file') > 0, 'No such file: %s', spm_mat);
    matlabbatch{e}.spm.stats.fmri_est.spmmat = cellstr(spm_mat);
    matlabbatch{e}.spm.stats.fmri_est.method.Classical = 1;


end

spm_jobman('run', matlabbatch);