function compute_contrasts(subs,exclude_subs,model,modelname)
%Compute secondlevel contrasts, depending on the GLM
%path business
if model == 1
    model_prefix = 'HRF';
elseif model == 2
    model_prefix = 'FIR';
end

% more path business
if strcmp(modelname,'split_onset_stimulus')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_stimulus')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_plateau')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_different_onsets_heat_pressure')
    modelstr = [model_prefix '_' modelname];
elseif strcmp(modelname,'complete_onset_stimulus_condition_contrast')
    modelstr = [model_prefix '_' modelname];
end

baseDir          = ['/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/']; % übergeordnet
FFXDir           = [baseDir 'derivatives/spm_firstlevel/' modelstr filesep];
RFXDir           = [baseDir 'derivatives/spm_secondlevel/'];
scriptDir        = ['/home/nold/Desktop/PEEP/fMRI/Analysis/fmri_analysis/second_level'];
%find mask, make sure it's there
mask     = fullfile(RFXDir, 'mask_secondlevel.nii');
%mask = fullfile(scriptDir,'neuromorphometrics.nii');


%deal with subjects and conditions
% exclude subjects
if ~isempty(exclude_subs)
    subs = subs(~ismember(subs,exclude_subs));
end

n_sub    = length(subs);
nSess       = 8;

if model == 1 && strcmp(modelname,'complete_onset_stimulus_condition_contrast')

    % define vectors for factors
    m_vec = [repmat([ones(1,6) -ones(1,6)],1,2)]; % modality
    i_vec = [repmat([-1 0 1],1,8)]; % intensity
    p_vec = [-ones(1,12) ones(1,12)]; %pharma
    e_vec = [repmat([-ones(1,3) ones(1,3)],1,4)]; % exercise

con_list = [
    {m_vec};    % heat>pressure
    {-m_vec};   % pressure>heat
    {i_vec};    % high>medum>low
    {-i_vec};
    {p_vec};    % NLX>SAL
    {-p_vec};
    {e_vec};    % high>low
    {-e_vec};
    {m_vec .* p_vec}; % modality*pharma
    {-(m_vec .* p_vec)};
    {m_vec .* i_vec}; % modality*intenisty
    {-(m_vec .* i_vec)};
    {m_vec .* e_vec}; % modality*exercise
    {-(m_vec .* e_vec)};
    {p_vec .* i_vec}; % pharma*intenisty
    {-(p_vec .* i_vec)};
    {p_vec .* e_vec}; % pharma*exercise
    {-(p_vec .* e_vec)};
    {i_vec .* e_vec}; % intensity*exercise
    {-(i_vec .* e_vec)};
    {m_vec .* p_vec .* e_vec}; % modality*pharma*exercise
    {-(m_vec .* p_vec .* e_vec)};
    {i_vec .* p_vec .* e_vec}; % intensity*pharma*exercise
    {-(i_vec .* p_vec .* e_vec)};
    {i_vec .* m_vec .* e_vec}; % intensity*modality*exercise
    {-(i_vec .* m_vec .* e_vec)};
    {i_vec .* m_vec .* p_vec}; % intensity*modality*pharma
    {-(i_vec .* m_vec .* p_vec)};
    
    ];
    

    %define names
    m = 'mod';
    i = 'int';
    e = 'ex';
    p = 'pharm';

    con_names = {
        m, 
        [m '_neg'],
        i, 
        [i '_neg'],
        p, 
        [p '_neg'],
        e, 
        [e '_neg'],

        ['int_' m '_' p],  
        ['int_' m '_' p '_neg'],
        ['int_' m '_' i],  
        ['int_' m '_' i '_neg'],
        ['int_' m '_' e],  
        ['int_' m '_' e '_neg'],
        ['int_' p '_' i],  
        ['int_' p '_' i '_neg'],
        ['int_' p '_' e],  
        ['int_' p '_' e '_neg'],
        ['int_' i '_' e],  
        ['int_' i '_' e '_neg'],

        ['int_' m '_' p '_' e],  
        ['int_' m '_' p '_' e '_neg'],
        ['int_' i '_' p '_' e],  
        ['int_' i '_' p '_' e '_neg'],
        ['int_' i '_' m '_' e],  
        ['int_' i '_' m '_' e '_neg'],
        ['int_' i '_' m '_' p],  
        ['int_' i '_' m '_' p '_neg'],

        };



    outdir = fullfile([RFXDir, model_prefix, '_', modelname]);
    spm_mat = fullfile(outdir, 'SPM.mat');
    assert(exist(spm_mat, 'file') > 0, 'No so SPM.mat: %s', spm_mat);


    matlabbatch = {};
    matlabbatch{1}.spm.stats.con.delete = 1;
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_mat);

    co = 1;

    matlabbatch{1}.spm.stats.con.consess{co}.fcon.name   = 'eoi';
    matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = eye(24);

    for count = 1:length(con_names)
        co = co + 1;
        matlabbatch{1}.spm.stats.con.consess{co}.tcon.name = con_names{count};
        matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec = con_list{count};
    end


% switch model
% 
%     case 1
%        % names = {'face_onsets', 'pred_fg', 'UCS'};
%     case 2
% 
%        modality = {'heat', 'pressure' };
%        intensity = {'30' , '50' , '70'};
%        %exercise = {'low', 'high'};
%        pharm_cond = {'nacl','nax'};
% 
% 
%        nBase = 12;
%        bins = 1:nBase;
%        s = 1;
%        names = {};
% 
%        % contrast names for each modality, intensity and bin
%        for p = 1:length(pharm_cond)
%            
%            for m = 1:length(modality)
% 
%                for i = 1:length(intensity)
% 
%                    for b = 1:length(bins)
% 
%                        names{s} = {[modality{m},intensity{i},'-',num2str(bins(b)),'_',pharm_cond{p}]};
%                        s = s + 1;
%                    end
% 
%                end
% 
%            end
%        end
% 
% 
% end
% 
% 
% 
% n_cons = length(names);
% con_joint = [repmat(eye(n_cons),1,nSess),ones(n_cons,length(subIDs))./length(subIDs)];
% 
% 
% sl_dir = [baseDir filesep 'derivatives' filesep 'spm_secondlevel'];
% 
% outdir = fullfile(sl_dir, modelname);
% spm_mat = fullfile(outdir, 'SPM.mat');
% assert(exist(spm_mat, 'file') > 0, 'No so SPM.mat: %s', spm_mat);
% 
% co = 0;
% matlabbatch = {};
% matlabbatch{1}.spm.stats.con.delete = 1;
% matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_mat);



% if model == 1
% 
%     %positive and negative tunings
%     co = co + 1;
%     matlabbatch{1}.spm.stats.con.consess{co}.tcon.name = '';
%     matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec = con_joint(2, :);
%     co = co + 1;
%     matlabbatch{1}.spm.stats.con.consess{co}.tcon.name = '';
%     matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec = - con_joint(2, :);

% elseif model == 2
% 
%         counter_0 = 0;
%         counter_1 = 0;
% 
%     
%    %t-contrast for each bin
%    %==========================================
% 
%     for idx = 1:n_cons
% 
%         co = co + 1;
%         matlabbatch{1}.spm.stats.con.consess{co}.tcon.name = names{idx}{:};
%         matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec = [zeros(1,counter_0) 1 zeros(1,n_cons - (counter_0 +1)) (ones(1,length(subIDs))*1./(length(subIDs)))] ;
% 
%         counter_0 = counter_0 + 1;
%         counter_1 = counter_1 + 1;
%     end

%     %F- contrast effect of interest
%     %==================================================
%     co = co + 1;
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.name   = 'eoi';
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = con_joint; 


%     %F-contrast for each modality and intensity
%     %==================================================
%  
%     clear names idx counter_1 counter_0 s
%     s = 1;
% 
%     for m = 1:length(modality)
%         names{s} = [modality{m}];
%         s = s + 1;
%     end
% 
%     counter_0 = 0;
%     counter_1 = 1;
% 
%     for idx = 1:length(names)
%         co = co + 1;
%         matlabbatch{1}.spm.stats.con.consess{co}.fcon.name = names{idx};
%         matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = [zeros(nBase,counter_0*nBase) eye(nBase) zeros(nBase,(n_cons - (counter_1*nBase))) ones(nBase,length(subIDs))./length(subIDs)];
%         counter_0 = counter_0 + 1;
%         counter_1 = counter_1 + 1;
%     end
%     
% 
%     %F-contrast for each modality and intensity
%     %==================================================
% 
%     clear names idx counter_1 counter_0 s
%     s = 1;
%  
%        for m = 1:length(modality)
%            for i = 1:length(intensity)
%                names{s} = [modality{m},'-',intensity{i}];
%                s = s + 1;
%            end
%        end
% 
%        counter_0 = 0;
%        counter_1 = 1; 
% 
%     for idx = 1:length(names)
%         co = co + 1;
%         matlabbatch{1}.spm.stats.con.consess{co}.fcon.name = names{idx};
%         matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = [zeros(nBase,counter_0*nBase) eye(nBase) zeros(nBase,(n_cons - (counter_1*nBase))) ones(nBase,length(subIDs))./length(subIDs)];
%         counter_0 = counter_0 + 1;
%         counter_1 = counter_1 + 1;
%     end
%     
% 
%     %F-contrast for pharm cond
%     %==================================================
%     co = co + 1;
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.name   = 'eoi_nax_nacl';
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = con_joint; 
%     
%      
%     %F-contrast for pharm cond and modality
%     %==================================================
%     co = co + 1;
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.name = 'heat_nax_nacl';
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec =[eye(n_cons/2),zeros(n_cons/2),ones(n_cons/2,length(subIDs))./length(subIDs)];
%     co = co + 1;
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.name = 'pressure_nax_nacl';
%     matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec =[zeros(n_cons/2),eye(n_cons/2),ones(n_cons/2,length(subIDs))./length(subIDs)];
%     %F-contrast for each modality and intensity and exercise intensity
 %=================================================================
 %   clear names idx counter_1 counter_0 s
 %   s = 1;
 
 %      for m = 1:length(modality)
%            for i = 1:length(intensity)
%                for e = 1:length(exercise)
%                    names{s} = [modality{m},'-',intensity{i},'-',exercise{e}];
%                    s = s + 1;
%                end
%            end
%        end
% 
%        counter_0 = 0;
%        counter_1 = 1; 
% 
%     for idx = 1:length(names)
%         co = co + 1;
%         matlabbatch{1}.spm.stats.con.consess{co}.fcon.name = names{idx};
%         matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = [zeros(nBase,counter_0*nBase) eye(nBase) zeros(nBase,(n_cons - (counter_1*nBase))) ones(nBase,length(subIDs))./length(subIDs)];
%         counter_0 = counter_0 + 1;
%         counter_1 = counter_1 + 1;
%     end
%     

% end


spm_jobman('run', matlabbatch);

end
