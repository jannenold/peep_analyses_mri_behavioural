function write_RFX_model(subs,exclude_subs,model,modelname,skern)
%function write_secondlevel_model(model)
%
%Specify the secondlevel model and write SPM.mat
%females_only = 0;
%males_only = 1;

%template for con images
file_tmp1= 's%dw_nlco_dartelcon_%04d.nii';
templateDir ='/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/spm_bids_pipeline/templates/';
estimate_contrasts_heat_pressure = 1;

%path business
if model == 1
    model_prefix = 'HRF';
elseif model == 2
    model_prefix = 'FIR';
end

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
    modelstr = [model_prefix '_' modelname '_6mm_heat_pressure'];
end

[path]  = get_study_specs;
FFXDir           = [path.firstlevelDir filesep modelstr];
RFXDir           = [path.secondlevelDir];

%find mask, make sure it's there
%mask     = fullfile(RFXDir, 'mask_secondlevel.nii');
mask = fullfile(templateDir,'neuromorphometrics.nii');

assert(exist(mask, 'file') > 0, 'No secondlevel mask.');

%deal with subjects and conditions
% exclude subjects

if ~isempty(exclude_subs)
    subs = subs(~ismember(subs,exclude_subs));
end

% if wanted: do only females/males
%fid = fopen ('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/participants.tsv');
%C = textscan (fid, '%s %d %d %d %d %d %f %f', 'HeaderLines', 1);
%fclose (fid);
%gender_vec = double(C{3});

%if females_only
%    subs = subs(gender_vec == 1);
%elseif males_only
%    subs = subs(gender_vec == -1);
%end


n_sub    = length(subs);

if model == 1 && strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm']) && ~estimate_contrasts_heat_pressure


%     names = {'pain_nacl','pain_nax','heat_nacl','heat_nax','pressure_nacl','pressure_nax',...
%         'pain_exHi_nacl','pain_exLo_nacl','pain_exHi_nax','pain_exLo_nax',...
%         'heat_exHi_nacl','heat_exLo_nacl','heat_exHi_nax','heat_exLo_nax',...
%         'pressure_exHi_nacl','pressure_exLo_nacl','pressure_exHi_nax', 'pressure_exLo_nax',...
%         'heat_nax>nacl','pressure_nax>nacl','heat_exHi>exLo_nacl','heat_exLo>exHi_nacl',...
%         'pressure_exHi>exLo_nacl','pressure_exLo>exHi_nac'};%,'heat>pressure_nacl','pressure>heat_nacl','rating'};

    % names = {'pain_nacl','pain_nax','heat_nacl','heat_nax'};%,'pressure_nacl','pressure_nax'};%,'heat_nacl_early','heat_nax_early','heat_nacl_late','heat_nax_late'};%'pain_nacl>nax','pain_nax>nacl','heat_nacl>nax','heat_nax>nacl','pressure_nacl>nax','pressure_nax>nacl', 'heat_nacl_int','heat_nax_int','pressure_nacl_int','pressure_nax_int','pressure>heat_nacl','pressure>heat_nax','heat>pressure_nacl','heat>pressure_nax','rating','pain_exHI>LO_nacl','pain_exHI>LO_nax','pain_exLO>HI_nacl','pain_exLO>HI_nax','heat_exHI>LO_nacl','heat_exHI>LO_nax','heat_exLO>HI_nacl','heat_exLO>HI_nax','pressure_exHI>LO_nacl','pressure_exHI>LO_nax','pressure_exLO>HI_nacl','pressure_exLO>HI_nax'};
    names = {...
        'pain_sal',...
        'pain_nlx',...
        'heat_sal',...
        'heat_nlx',...
        'pressure_sal',...
        'pressure_nlx',...
        'heat_nlx>heat_sal',...
        'heat_sal>heat_nlx',...
        'pressure_nlx>pressure_sal',...
        'pressure_sal>pressure_nlx',...
        'heat_sal>pressure_sal',...
        'pressure_sal>heat_sal',...
        'heat_nlx>pressure_nlx',...
        'pressure_nlx>heat_nlx',...
        'heat_int_sal',...
        'heat_int_nlx',...
        'pressure_int_sal',...
        'pressure_int_nlx',...
        'heat_exLo_sal',...
        'heat30_exLo_sal',...
        'heat50_exLo_sal',...
        'heat70_exLo_sal',...
        'heat_exHi_sal',...
        'heat30_exHi_sal',...
        'heat50_exHi_sal',...
        'heat70_exHi_sal',...
        'heat_exLo_nlx',...
        'heat30_exLo_nlx',...
        'heat50_exLo_nlx',...
        'heat70_exLo_nlx',...
        'heat_exHi_nlx',...
        'heat30_exHi_nlx',...
        'heat50_exHi_nlx',...
        'heat70_exHi_nlx',...
        'pressure_exLo_sal',...
        'pressure30_exLo_sal',...
        'pressure50_exLo_sal',...
        'pressure70_exLo_sal',...
        'pressure_exHi_sal',...
        'pressure_exLo_nlx',...
        'pressure30_exLo_nlx',...
        'pressure50_exLo_nlx',...
        'pressure70_exLo_nlx',...
        'pressure_exHi_nlx',...
        'heat_exHi>exLo_sal',...
        'heat_exLo>exHi_sal',...
        'pressure_exHi>exLo_sal',...
        'pressure_exLo>exHi_sal',...
        'heat_exHi>exLo_nlx',...
        'heat_exLo>exHi_nlx',...
        'pressure_exHi>exLo_nlx',...
        'pressure_exLo>exHi_nlx',...
        'heat70_exHi>exLo_sal',...
        'heat70_exLo>exHi_sal',...
        'heat70_exHi>exLo_nlx',...
        'heat70_exLo>exHi_nlx',...
        'heat50_exHi>exLo_sal',...
        'heat50_exLo>exHi_sal',...
        'heat50_exHi>exLo_nlx',...
        'heat50_exLo>exHi_nlx',...
        'heat30_exHi>exLo_sal',...
        'heat30_exLo>exHi_sal',...
        'heat30_exHi>exLo_nlx',...
        'heat30_exLo>exHi_nlx',...
        'pressure70_exHi>exLo_sal',...
        'pressure70_exLo>exHi_sal',...
        'pressure70_exHi>exLo_nlx',...
        'pressure70_exLo>exHi_nlx',...
        'pressure50_exHi>exLo_sal',...
        'pressure50_exLo>exHi_sal',...
        'pressure50_exHi>exLo_nlx',...
        'pressure50_exLo>exHi_nlx',...
        'pressure30_exHi>exLo_sal',...
        'pressure30_exLo>exHi_sal',...
        'pressure30_exHi>exLo_nlx',...
        'pressure30_exLo>exHi_nlx',...
        'inter_modality_treatment',...
        'inter_intensity_treatment',...
        'inter_intensity_treatment_heat',...
        'inter_intensity_treatment_pressure',...
        'inter_intensity_modality_sal',...
        'inter_intensity_modality_nlx',...
        'heat30_sal',...
        'heat50_sal',...
        'heat70_sal',...
        'pressure30_sal',...
        'pressure50_sal',...
        'pressure70_sal',...
        'heat30_nlx',...
        'heat50_nlx',...
        'heat70_nlx',...
        'pressure30_nlx',...
        'pressure50_nlx',...
        'pressure70_nlx',...
        'heat',...
        'pressure',...
        'inter_exercise_treatment_heat',...
        'inter_exercise_treatment_pressure',...
        'inter_exercise_treatment_heat70',...
        'rating'};

    n_cons = length(names);
    co_c = 0;
    
    %collect con images over all subjects and conditions
    for e = 1:n_cons

        outdir = fullfile(RFXDir, modelstr,names{e});

        if ~exist(outdir)
            mkdir(outdir);
        end

        matlabbatch{e}.spm.stats.factorial_design.dir = cellstr(outdir);

        co_c = co_c + 1;
        all_niftis = cell(length(n_sub), 1);


        for s = 1:n_sub

            path = fullfile(FFXDir, sprintf('sub-%02d',subs(s)));
            templ = sprintf(file_tmp1, skern, e);
            nifti = spm_select('ExtFPList', path, templ);
            assert(size(nifti, 1) == 1);
            all_niftis{s,1} = cellstr(nifti);

        end

        % load covariate (FTP value z-standardised)
        fid = fopen ('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/participants.tsv'); 
        C = textscan (fid, '%s %d %d %d %d %d %f %f', 'HeaderLines', 1); 
        fclose (fid);
        ftp_vec = double(C{5});
        pwc_vec = C{8};
        expect_pain_rating = C{7};
        gender_vec = double(C{3});


       
        all_niftis = {cat(1,all_niftis{:})};
        matlabbatch{e}.spm.stats.factorial_design.des.t1.scans = all_niftis{:};
        %matlabbatch{e}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.cov.c = [gender_vec];
        matlabbatch{e}.spm.stats.factorial_design.cov.cname = 'sex';
        matlabbatch{e}.spm.stats.factorial_design.cov.iCFI = 1;
        matlabbatch{e}.spm.stats.factorial_design.cov.iCC = 1;
        matlabbatch{e}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.em = {mask};
        matlabbatch{e}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.glonorm = 1;
    end

elseif model == 1 && strcmp(modelname,['complete_onset_stimulus_' num2str(skern) 'mm']) && estimate_contrasts_heat_pressure

    names = {...
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
        'pressure70_exLo_sal',...
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

    n_cons = length(names);
    co_c = 0;
    
    %collect con images over all subjects and conditions
    for e = 1:n_cons

        outdir = fullfile(RFXDir, modelstr,names{e});

        if ~exist(outdir)
            mkdir(outdir);
        end

        matlabbatch{e}.spm.stats.factorial_design.dir = cellstr(outdir);

        co_c = co_c + 1;
        all_niftis = cell(length(n_sub), 1);


        for s = 1:n_sub

            path = fullfile(FFXDir, sprintf('sub-%02d',subs(s)));
            templ = sprintf(file_tmp1, skern, e);
            nifti = spm_select('ExtFPList', path, templ);
            assert(size(nifti, 1) == 1);
            all_niftis{s,1} = cellstr(nifti);

        end

       
        all_niftis = {cat(1,all_niftis{:})};
        matlabbatch{e}.spm.stats.factorial_design.des.t1.scans = all_niftis{:};
        matlabbatch{e}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.em = {mask};
        matlabbatch{e}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.glonorm = 1;
    end

elseif model == 1 && strcmp(modelname,['split_onset_stimulus']) && estimate_contrasts_heat_pressure

    names = {'heat_exlo_sal_early',...
        'heat_exlo_sal_late',...
        'pressure_exlo_sal_early',...
        'pressure_exlo_sal_late',....
        'heat_exlo_sal_early_late',...
        'pressure_exlo_sal_early_late',...
        'early_exlo_sal_heat_pressure',...
        'late_exlo_sal_heat_pressure'};

    n_cons = length(names);
    co_c = 0;
    
    %collect con images over all subjects and conditions
    for e = 1:n_cons

        outdir = fullfile(RFXDir, modelstr,names{e});

        if ~exist(outdir)
            mkdir(outdir);
        end

        matlabbatch{e}.spm.stats.factorial_design.dir = cellstr(outdir);

        co_c = co_c + 1;
        all_niftis = cell(length(n_sub), 1);


        for s = 1:n_sub

            path = fullfile(FFXDir, sprintf('sub-%02d',subs(s)));
            templ = sprintf(file_tmp1, skern, e);
            nifti = spm_select('ExtFPList', path, templ);
            assert(size(nifti, 1) == 1);
            all_niftis{s,1} = cellstr(nifti);

        end

       
        all_niftis = {cat(1,all_niftis{:})};
        matlabbatch{e}.spm.stats.factorial_design.des.t1.scans = all_niftis{:};
        matlabbatch{e}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.em = {mask};
        matlabbatch{e}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.glonorm = 1;
    end

elseif  model == 1 && strcmp(modelname,'split_onset_stimulus') && ~estimate_contrasts_heat_pressure


    names = {'heat_nacl_early','heat_nacl_late','pressure_nacl_early','pressure_nacl_late','heat_nlx_early','heat_nlx_late','pressure_nlx_early','pressure_nlx_late','heat_early_nlx_sal',...
        'heat_early_sal_nlx','heat_late_nlx_sal','heat_late_sal_nlx','pressure_early_nlx_sal', 'pressure_early_sal_nlx','pressure_late_nlx_sal','pressure_late_sal_nlx','heat_int_sal_early',...
        'heat_int_sal_late','heat_int_nlx_early','heat_int_nlx_late','pressure_int_sal_early','pressure_int_sal_late','pressure_int_nlx_early','pressure_int_nlx_late',...
        'inter_treatment_early_late_heat','inter_treatment_early_late_pressure','inter_intensity_treatment_heat_early','inter_intensity_treatment_heat_late',...
        'inter_intensity_treatment_pressure_early','inter_intensity_treatment_pressure_late','inter_modality_treatment_early','inter_modality_treatment_late',...
        'heat_exLo>exHi_sal_early','heat_exLo>exHi_sal_late','heat_exLo>exHi_nlx_early','heat_exLo>exHi_nlx_late'};
    n_cons = length(names);
    co_c =0;
    %collect con images over all subjects and conditions
    for e = 1:n_cons

        outdir = fullfile(RFXDir, modelstr,names{e});

        if ~exist(outdir)
            mkdir(outdir);
        end

        matlabbatch{e}.spm.stats.factorial_design.dir = cellstr(outdir);

        co_c = co_c + 1;
        all_niftis = cell(length(n_sub), 1);


        for s = 1:n_sub

            path = fullfile(FFXDir, sprintf('sub-%02d',subs(s)));
            templ = sprintf(file_tmp1, skern, e);
            nifti = spm_select('ExtFPList', path, templ);
            assert(size(nifti, 1) == 1);
            all_niftis{s,1} = cellstr(nifti);

        end

        fid = fopen ('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/participants.tsv');
        C = textscan (fid, '%s %d %s %d %d %d %f %f', 'HeaderLines', 1);
        fclose (fid);
        ftp_vec = double(C{5});
        pwc_vec = C{8};
        expect_pain_rating = C{7};
        gender = char(C{3});

        t = strrep(gender.','f','1');
        t = strrep(t,'m','0');
        gender_vec = str2num(t.');


        all_niftis = {cat(1,all_niftis{:})};
        matlabbatch{e}.spm.stats.factorial_design.des.t1.scans = all_niftis{:};
        matlabbatch{e}.spm.stats.factorial_design.cov.c = [pwc_vec];
        matlabbatch{e}.spm.stats.factorial_design.cov.cname = 'rftp';
        matlabbatch{e}.spm.stats.factorial_design.cov.iCFI = 1;
        matlabbatch{e}.spm.stats.factorial_design.cov.iCC = 1;
        matlabbatch{e}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{e}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{e}.spm.stats.factorial_design.masking.em = {mask};
        matlabbatch{e}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{e}.spm.stats.factorial_design.globalm.glonorm = 1;
    end

elseif model == 1 && strcmp(modelname,'complete_onset_stimulus_condition_contrast') %% Flexible Factorial


    outdir = fullfile(RFXDir, modelstr);
    if ~exist(outdir)
        mkdir(outdir);
    end


    matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'GROUP';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'SUBJECT';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'CONDITION';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;

    %GROUP SUBJECT DAY CONDITION
    allfiles  = [];
    imat      = [];
    inc_su = 1;
    groups = {1};
    cond_use = 1:24;
    warp_con = 1;
    file_templ   = 's%dw_nlco_dartelcon_%04d.nii';


    for gr = 1:size(groups,2)
        for su = 1:size(subs,2)
            for co = 1:size(cond_use,2)
                name        = sprintf('sub-%02d',subs(su));
                %a_dir       = [base_dir name filesep anadirname filesep];
                %s_string    = sprintf('s%s',sm_str);
                %if skern == 0;s_string = '';end
                %if warp_beta
                %    swcon_file = [a_dir sprintf(['%s%s' beta_temp], s_string, n_type,simple_beta(co))];
                %end
                if warp_con
                    swcon_file = fullfile(FFXDir, name, sprintf(file_templ,skern,co));
                end
                allfiles    = strvcat(allfiles,swcon_file);
                mat_entry   = [1 gr subs(su) co];
                imat        = [imat; mat_entry];
            end

        end
    end

    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans    = cellstr(allfiles);
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix  = imat;

    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum  = [3];
    %matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum  = [2];
    %matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [1 3];


    matlabbatch{1}.spm.stats.factorial_design.cov                  = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov            = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none   = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im           = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em           = {mask};


   
elseif model == 1 && strcmp(modelname,'split_onset_stimulus_fact') && estimate_contrasts_heat_pressure %% Flexible Factorial

FFXDir           = [path.firstlevelDir filesep 'HRF_split_onset_stimulus_6mm_heat_pressure'];
    outdir = fullfile(RFXDir, ['FACT_',modelstr]);
    if ~exist(outdir)
        mkdir(outdir);
    end


    matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'GROUP';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'SUBJECT';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'CONDITION';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;

    %GROUP SUBJECT DAY CONDITION
    allfiles  = [];
    imat      = [];
    inc_su = 1;
    groups = {1};
    cond_use = 1:8;
    warp_con = 1;
    file_templ   = 's%dw_nlco_dartelcon_%04d.nii';


    for gr = 1:size(groups,2)
        for su = 1:size(subs,2)
            for co = 1:size(cond_use,2)
                name        = sprintf('sub-%02d',subs(su));
                %a_dir       = [base_dir name filesep anadirname filesep];
                %s_string    = sprintf('s%s',sm_str);
                %if skern == 0;s_string = '';end
                %if warp_beta
                %    swcon_file = [a_dir sprintf(['%s%s' beta_temp], s_string, n_type,simple_beta(co))];
                %end
                if warp_con
                    swcon_file = fullfile(FFXDir, name, sprintf(file_templ,skern,co));
                end
                allfiles    = strvcat(allfiles,swcon_file);
                mat_entry   = [1 gr subs(su) co];
                imat        = [imat; mat_entry];
            end

        end
    end

    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans    = cellstr(allfiles);
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix  = imat;

    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum  = [3];
    %matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum  = [2];
    %matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [1 3];


    matlabbatch{1}.spm.stats.factorial_design.cov                  = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov            = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none   = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im           = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em           = {mask};



elseif model == 2

    nConds = 6;
    nBase = 12;
    co_c = 0; %contrast counter

    %collect con images over all subjects and conditions
    e = 1;

    N_conimages = nBase*nConds;
    nco = 1;
    outdir = fullfile(RFXDir, modelstr);

    if ~exist(outdir)
        mkdir(outdir);
    end


    co_c = co_c + 1;
    niftis = cell(N_conimages, 1);

    for s = 1:n_sub

        for cons = 1:N_conimages
            path = fullfile(FFXDir, sprintf('sub-%02d',subs(s)));
            temp2 = sprintf(file_tmp1, skern,cons);
            nifti = spm_select('FPList', path, temp2);
            assert(size(nifti, 1) == 1);
            niftis{cons} = nifti;
        end

        matlabbatch{e}.spm.stats.factorial_design.des.anovaw.fsubject(s).scans = niftis;
        matlabbatch{e}.spm.stats.factorial_design.des.anovaw.fsubject(s).conds = 1:N_conimages;


    end


    matlabbatch{e}.spm.stats.factorial_design.des.anovaw.dept  = 0;
    matlabbatch{e}.spm.stats.factorial_design.des.anovaw.variance  = 0;
    matlabbatch{e}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
    matlabbatch{e}.spm.stats.factorial_design.des.anovaw.ancova = 0;
    matlabbatch{e}.spm.stats.factorial_design.dir = cellstr(outdir);
    matlabbatch{e}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{e}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{e}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{e}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{e}.spm.stats.factorial_design.masking.em = {mask};
    matlabbatch{e}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{e}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{e}.spm.stats.factorial_design.globalm.glonorm = 1;




end




% run matlabbatch
spm_jobman('run',matlabbatch);







