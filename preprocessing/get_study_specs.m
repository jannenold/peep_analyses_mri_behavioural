function [path,vars,analysis,import] = get_study_specs


%% import definitions
import.prisma        = {{24170, 24192},{24179, 24198},{24194, 24217},{24211, 24229},{24280, 24302},{24312, 24336},{24324, 24350},{24335, 24355},{24368, 24390},...
                       {24365, 24449},{24414, 24422},{24521, 24539},{24507, 24525},{24453, 24468},{24470, 24572},{24514, 24566},{24528, 24544},{24547, 24570},...
                       {24542, 24686},{24557, 24571},{24625, 24627},{24573, 24587},{24578, 24613},{24609, 24622},{24642, 24652},{24633, 24648},{24644, 24675},...
                       {24646, 24677},{24697, 24725},{24680, 24690},{24699, 24733},{24720, 24727},{24772, 24776},{24747, 24765},{24800, 24826},{24778, 24791},...
                       {24795, 24804},{24808, 24819},{24811, 24823},{24824, 24828}};
import.prisma_no     = [1,              2,              3,              4,          6,              7,              8,              9,          11,...
                        13,             14,             15,             17,         18,             19,             21,             22,         23,...
                        24,             25,             26,             27,         28,             29,             30,             31,         32,...
                        34,             35,             36,             38,         39,             40,             41,             43,         44,...
                        45,             46,             47,             48];

% import prisma sub-26: anatomical scan on another day 
%import.prisma = {{24606}};
%import.prisma_no = [26];
                    
import.user          = 'nold';
import.server        = 'revelations.nin.uke.uni-hamburg.de';

import.data(1).dir        = 'func'; 
import.data(1).type       = 'epi';
import.data(1).seq        = 'ep2d_bold, 2.0mm3, mb2, fMRI  '; %protocol name (trailing space makes it unique) 
import.data(1).cond       = 'n > 300'; % heuristic to get only valid runs (e.g. more than 1000 volumes)

import.data(2).dir        = 'anat'; % valid BIDS dir name
import.data(2).type       = 'T1w'; % valid BIDS file name
import.data(2).seq        = 'ninFLASH_v14A_df, mprage, defa-SAT_DEFA ';
import.data(2).cond       = 'n == 240'; % heuristic to get only valid runs (e.g. exactly 240 slices)
   
import.dummies            = 4; 

%% path definitions

path.baseDir     = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN';
path.templateDir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/spm_bids_pipeline/templates/';

path.rawDir          = fullfile(path.baseDir, 'rawdata'); % that should be the folder with the zipped not-to-touch data sets
path.derivDir        = fullfile(path.baseDir, 'derivatives');
path.preprocDir      = fullfile(path.baseDir, 'derivatives', 'spm_preprocessing');
path.firstlevelDir   = fullfile(path.baseDir, 'derivatives', 'spm_firstlevel');
path.secondlevelDir  = fullfile(path.baseDir, 'derivatives', 'spm_secondlevel');

%% vars definitions

% various predefined names (do not change)
vars.skullStripID    = 'skull-strip-T1.nii';
vars.T1maskID        = 'brain_mask.nii';
vars.templateID      = 'cb_Template_%d_Dartel.nii';
vars.templateT1ID    = 'cb_Template_T1.nii';
%vars.groupMaskID     = 'brainmask.nii';
vars.groupMaskID     = 'neuromorphometrics.nii';
%% this need to be adapted to your study / computer--------------
vars.max_procs   = 8;
vars.task        = 'peep';
vars.nRuns       = 4;
vars.nSess       = 2;
% get info for slice timing correction
vars.sliceTiming.so    = [1720, 1660, 1603, 1543, 1483, 1423, 1365, 1305, 1245, 1188, 1128, 1068, 1008, 950, 890, 830, 770, 713, 653, 593, 535, 475, 415, 355, 298, 238, 178, 118, 60,  0,...
                          1720, 1660, 1603, 1543, 1483, 1423, 1365, 1305, 1245, 1188, 1128, 1068, 1008, 950, 890, 830, 770, 713, 653, 593, 535, 475, 415, 355, 298, 238, 178, 118, 60,  0]; % in ms
vars.sliceTiming.tr       = 1.8; % in s
vars.sliceTiming.nslices  = 60;
vars.sliceTiming.refslice = 900;



%% analysis definitions

% little routine to generate cond_names as in the TSV files
treatment = {'sal','nlx'};
exercise  = {'lo','hi'};
modality  = {'heat','pressure'};
intensity = {'30','50','70'};

cnt = 1;
for tt = 1:numel(treatment)
    for ee = 1:numel(exercise)
        for mm = 1:numel(modality)
            for ii = 1:numel(intensity)
                cond_names{cnt} = [treatment{tt} '_' exercise{ee} '_' modality{mm} '_' intensity{ii}];
                cnt = cnt + 1;
            end
        end
    end
end



analysis.all_subs  = [1:4 6:9 11 13:15 17:19 21:32 34:36 38:41 43:48]; % very all;
sporty             = [0 1 0 1 0 0 0 1 0 0 1 1 0 1 1 1 1 0 1 0 0 0 1 1 1 0 0 1 0 1 0 0 1 1 1 1 0 1 0 0];
expect             = [1 1 0 1 0 0 1 1 1 1 0 1 0 1 0 0 0 1 1 0 1 0 1 0 1 0 1 0 1 1 1 0 1 0 0 1 1 0 0 0];
single_group       = ones(size(expect));

% analysis.group_ind  = single_group; % 1 for everybody ie 1 group for one-sample t-test
% analysis.group_weights = [1];
% analysis.group_names   = {'All'};

analysis.group_ind  = sporty+1; %index 1 and 2
analysis.group_weights = [1 1;1 -1;-1 1];
analysis.group_names   = {'All','Sporty>Lazy','Lazy>Sporty'};



analysis.parallel          = 0;
analysis.noise_corr        = ['mov24_physio'];
%analysis.noise_corr        = ['physio'];
analysis.cvi               = 'none'; % 'AR(1)'  'FAST' 'none'
analysis.shift             = 0; % shift all onsets
analysis.skernel           = 6; % smoothing kernel
analysis.wls               = 0;
analysis.bs                = 1;
analysis.concatenate       = 0;
analysis.cond_names        = cond_names; % they need to be EXACTLY the same as in the TSV files!!!
%what to do
analysis.do_model          = 0;
analysis.do_est            = 0;
analysis.do_vasa           = 0;
analysis.do_cons           = 0;
analysis.do_correct_vasa   = 0;
analysis.do_warp           = 0;
analysis.do_smooth         = 0;

analysis.do_fact           = 1;
analysis.do_fact_con       = 1;

analysis.do_one_t          = 0;


do_hrf = 1;
do_fir = 0;

if do_hrf
    analysis.max_procs         = 12;
    
    [analysis.t_con, analysis.t_con_names] = get_hrf_cons;
    
    analysis.ana               = 2; %hrf
    analysis.n_base            = 1;
    analysis.name_ana          = 'peep_hrf';
    
elseif do_fir
    analysis.max_procs         = 8; % big design matrix ...
    analysis.t_con             = [];
    analysis.t_con_names       = {};
    analysis.ana               = 1; %fir
    analysis.n_base            = 12;
    analysis.name_ana          = 'peep_fir';
end
end

function [t_con, t_con_names] = get_hrf_cons(cond_names)
t_vec     = [-ones(1,12) ones(1,12)]; % sal nlx
e_vec     = repmat([-ones(1,6) ones(1,6)],1,2); %low hi
m_vec     = repmat([-ones(1,3) ones(1,3)],1,4); % heat pressure
i_vec     = repmat([-1 0 1],1,8); % 30 50 70

% simple main effects

t_con_names = {'inter_int_nlx_heat70','neg_inter_int_nlx_heat70','int_heat','neg_int_heat','int_press','neg_int_press','inter_int_mod','neg_inter_int_mod',...
               'inter_int_nlx','neg_inter_int_nlx','inter_int_nlx_heat','neg_inter_int_nlx_heat',...
               'inter_int_nlx_press','neg_inter_int_nlx_press','inter_exc_heat','neg_inter_exc_heat','inter_exc_press','neg_inter_exc_press',...
               'inter_exc_heat_treat','neg_inter_exc_heat_treat','inter_exc_heat70_treat','neg_inter_exc_heat70_treat','inter_exc_press70_treat','neg_inter_exc_press70_treat'};

t_con = [
         [(m_vec<0) .* (i_vec>0) .* t_vec];...
         -[(m_vec<0) .* (i_vec>0) .* t_vec];...
         [(m_vec<0) .* i_vec];...
         -[(m_vec<0) .* i_vec];...
         [(m_vec>0) .* i_vec];...
         -[(m_vec>0) .* i_vec];...
         [m_vec .* i_vec];...
         -[m_vec .* i_vec];...
         [t_vec .* i_vec];...
         -[t_vec .* i_vec];...
         [(m_vec<0) .*t_vec .* i_vec];...
         -[(m_vec<0) .*t_vec .* i_vec];...
         [(m_vec>0) .*t_vec .* i_vec];...
         -[(m_vec>0) .*t_vec .* i_vec];...         
         [(m_vec<0)  .* i_vec .* e_vec];...
         -[(m_vec<0)  .* i_vec .* e_vec];...        
         [(m_vec>0)  .* i_vec .* e_vec];...
         -[(m_vec>0)  .* i_vec .* e_vec];...        
         [(m_vec<0)  .* i_vec .* e_vec.* t_vec];...
         -[(m_vec<0)  .* i_vec .* e_vec.* t_vec];...        
         [(m_vec<0)  .* (i_vec>0) .* e_vec.* t_vec];...
         -[(m_vec<0)  .* (i_vec>0) .* e_vec.* t_vec];...        
         [(m_vec>0)  .* (i_vec>0) .* e_vec.* t_vec];...
         -[(m_vec>0)  .* (i_vec>0) .* e_vec.* t_vec];...        
         
         ];

end