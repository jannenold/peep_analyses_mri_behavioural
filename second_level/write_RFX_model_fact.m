function write_RFX_model_fact
%function write_secondlevel_model(model)
%
%Specify the secondlevel model and write SPM.mat
do_con = 1;
clear all;
subs          = [1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48];
exclude_subs    = [1,4,30:32,46];

if ~isempty(exclude_subs)
    subs = subs(~ismember(subs,exclude_subs));
end
matlabbatch = [];
allFiles = [];
assembledCons = [];
condUse = [1:30];
%condUse = [1:12];
filesByCon = [];
file_tmp1= 's%dw_nlco_dartelcon_%04d.nii';
%FFXDir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_firstlevel/HRF_split_onset_stimulus_6mm_heat_pressure';
FFXDir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_firstlevel/HRF_complete_onset_stimulus_6mm_heat_pressure';
skern = 6;
out_dir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_complete_onset_not_all_subs';
%out_dir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/FACT_heat_pressure_split_onset';
templateDir ='/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/spm_bids_pipeline/templates/';
mask = fullfile(templateDir,'neuromorphometrics.nii');

for co = 1:size(condUse,2)
    filesByCon{co} = [];
end

for g = 1:numel(subs)
    name       = sprintf('sub-%02d',subs(g));
    allFiles = [];
    assembledCons = [];

    for co = 1:size(condUse,2)

        swbeta_templ      = sprintf(file_tmp1, skern, condUse(co)); % TO USE CONTRASTS
        allFiles = strvcat(allFiles,[FFXDir filesep name filesep swbeta_templ]); % all cons or betas included
        filesByCon{co} = strvcat(filesByCon{co},[FFXDir filesep name filesep swbeta_templ]);
        assembledCons = [assembledCons condUse(co)];
    end


end


%% --------------------- MODEL SPECIFICATION --------------------- %%

for co = 1:size(condUse,2)
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(co).scans = cellstr(filesByCon{co}); % presumably
end
matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {mask};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% --------------------- MODEL ESTIMATION --------------------- %%

matlabbatch{2}.spm.stats.fmri_est.spmmat = {[out_dir '/SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch);

%% --------------------- CONTRASTS --------------------- %%

matlabbatch = [];

clear SPM; load([out_dir '/SPM.mat']); %should exist by now
matlabbatch{1}.spm.stats.con.spmmat = {[out_dir '\SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 1;

co = 1;
% 
% t_con_names = {'heat_exlo_sal_early',...% heat early
%     'heat_exlo_sal_late',...            % heat late
%     'pressure_exlo_sal_early',...       % pressure early
%     'pressure_exlo_sal_late',....       % pressure late
%     'heat_exlo_sal_early_late',...      % heat early > late
%     'pressure_exlo_sal_early_late',...  % pressure early > late
%     'early_exlo_sal_heat_pressure',...  % early pain heat > pressure
%     'late_exlo_sal_heat_pressure',...     % later pain heat > pressure
%     'heat70_exlo_sal_early_late',...      % heat early > late
%     'pressure70_exlo_sal_early_late',...  % pressure early > late
%     'early_exlo_sal_heat70_pressure70',...  % early pain heat > pressure
%     'late_exlo_sal_heat70_pressure70'};


       t_con_names = {...
            'heat_exLo_sal',...
            'heat30_exLo_sal',...
            'heat50_exLo_sal',...
            'heat70_exLo_sal',...
            'heat_exLo_nlx',...
            'heat30_exLo_nlx',...
            'heat50_exLo_nlx',...
            'heat70_exLo_nlx',...
            'heat_nlx>sal_exLo',...
            'heat70_nlx>sal_exLo',...
            'heat_int_sal_exLo',...
            'inter_int_treat_heat_exLo',...
            'pressure_exLo_sal',...
            'pressure30_exLo_sal',...
            'pressure50_exLo_sal',...
            'pressu1re70_exLo_sal',...
            'pressure_exLo_nlx',...
            'pressure30_exLo_nlx',...
            'pressure50_exLo_nlx',...
            'pressure70_exLo_nlx',...
            'pressure_nlx>sal_exLo',...
            'pressure70_nlx>sal_exLo',...
            'pressure_int_sal_exLo',...
            'inter_int_treat_pressure_exLo',...
            'inter_modality_treat_exLo',...
            'inter_modality_treat_70_exLo',...
            'heat_nlx>pressure_nlx_exLo',...
            'heat70_nlx>pressure70_nlx_exLo',...
            'inter_int_modality_sal_exLo',...
            'heat_sal>pressure_sal_exLo'};


t_con = eye(size(t_con_names,2));



% obtain t contrasts
for gco = 1:size(t_con_names,2)

    matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = t_con_names{gco};
    matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = t_con(gco,:);
    matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
    co = co + 1; % increment by 1
end



spm_jobman('run',matlabbatch);
