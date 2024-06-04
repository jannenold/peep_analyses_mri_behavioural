clear all; clc
subs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 
exclude_subs    = [8];
model           = 1; % 1 = HRF, 2 = FIR;
skern           = 6; % smoothing kernel from FFX
modelname       = ['complete_onset_stimulus_' num2str(skern) 'mm']; %full length, onset start stimulus

file_tmp1= 's%dw_nlco_dartelcon_0078.nii';



%path business
if model == 1
    model_prefix = 'HRF';
elseif model == 2
    model_prefix = 'FIR';
end

% more path business
if strcmp(modelname,'split_onset_stimulus')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm'])
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_plateau')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_different_onsets_heat_pressure')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_stimulus_condition_contrast')
    modelstr = [model_prefix '_' modelname];
end


baseDir          = ['/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/']; % übergeordnet
FFXDir           = [baseDir 'derivatives/spm_firstlevel/' modelstr filesep];
RFXDir           = [baseDir 'derivatives/spm_secondlevel/'];

%find mask, make sure it's there
mask     = fullfile(RFXDir, 'mask_secondlevel.nii');
%mask = fullfile(scriptDir,'neuromorphometrics.nii');

assert(exist(mask, 'file') > 0, 'No secondlevel mask.');

%deal with subjects and conditions
% exclude subjects

if ~isempty(exclude_subs)
    subs = subs(~ismember(subs,exclude_subs));
end

% load covariate (FTP value z-standardised)
fid = fopen ('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/participants.tsv');
C = textscan (fid, '%s %d %d %d %d %d %f %f', 'HeaderLines', 1);
fclose (fid);
pwc_vec = C{8};
gender_vec = double(C{3});



pwc_vec_male = pwc_vec(gender_vec==-1);
pwc_vec_new = [pwc_vec_male];


subs_male = subs(gender_vec==-1);
n_sub_male    = length(subs_male);


names = {'contrast'};

n_cons = length(names);
co_c = 0;

%collect con images over all subjects and conditions
for e = 1:n_cons

    outdir = fullfile(RFXDir, '1sampleTtest');

    if ~exist(outdir)
        mkdir(outdir);
    end

    matlabbatch{e}.spm.stats.factorial_design.dir = cellstr(outdir);

    co_c = co_c + 1;
    all_niftis_male = cell(length(subs_male), 1);



    for b = 1:n_sub_male

        path = fullfile(FFXDir, sprintf('sub-%02d',subs_male(b)));
        templ = sprintf(file_tmp1, skern);
        nifti = spm_select('ExtFPList', path, templ);
        assert(size(nifti, 1) == 1);
        all_niftis_male{b,1} = cellstr(nifti);

    end
    all_niftis_male = {cat(1,all_niftis_male{:})};


    matlabbatch{e}.spm.stats.factorial_design.des.t1.scans = all_niftis_male{:};
    matlabbatch{e}.spm.stats.factorial_design.cov.c = [gender_vew];
    matlabbatch{e}.spm.stats.factorial_design.cov.cname = 'ftp';
    matlabbatch{e}.spm.stats.factorial_design.cov.iCFI = 1;
    matlabbatch{e}.spm.stats.factorial_design.cov.iCC = 1;
    matlabbatch{e}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{e}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{e}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{e}.spm.stats.factorial_design.masking.em = {mask};
    matlabbatch{e}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{e}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{e}.spm.stats.factorial_design.globalm.glonorm = 1;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


end

% run matlabbatch
spm_jobman('run',matlabbatch);