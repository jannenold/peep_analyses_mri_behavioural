function BS_FFX_peep(allSubs,exclude,session,model,modelname)
%firstlevel pipeline based on preprocessing pipeline
%pdev_preprocessing_fmri_unix.m
%u.bromberg@uke.de 20211007
% adapted by janne nold 28.06.2022

addpath(genpath('/home/nold/Desktop/PEEP/fMRI/'));

if ~isempty(exclude)
    allSubs = allSubs(~ismember(allSubs,exclude));
end

%% call hostname and define paths and number of processing units
hostName    = char(getHostName(java.net.InetAddress.getLocalHost));
studyName   = 'peep';

if model == 1
    basisF          = 'HRF'; %'HRF'         %one of HRF | FIR
elseif model == 2
    basisF          = 'FIR'; %'HRF'         %one of HRF | FIR
end

switch hostName
    case 'revelations.nin.uke.uni-hamburg.de'
        baseDir          = ['/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/']; % übergeordnet
        rawDir            = [baseDir 'rawdata/'];
        derivDir         = [baseDir 'derivatives/spm_preprocessing/']; % untergeordnet
        FFXDir           = [baseDir 'derivatives/spm_firstlevel/' 'BS_' basisF '_' modelname filesep];
        nProc            = 8; % 8 parallel processes on revelations
        subDir            = {'anat','func','fmap'};

    otherwise
        error('host not recognized');
end

%% define general studyspecific values / names

%subjects / runs / conditions
allRuns     = 1:4;
nSess       = 8;
nPmod       = 0;%2;                %number parametric modulators
imageDim    = '4D';             %dimension of images to be processed (3D or 4D)

%first-level specific values
orthog          = 1;                    %as input for conditions (0 or 1)
skern           = 4;                    %for smooting
skernel         = repmat(skern,1,3);    %for smoothing
units           = 'scans';              %onsets are specified in 'scans' or 'secs'
cvi             = 'none';               % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1)'|'AR(1) + w'
hpf             = 128;% 360;           % low frequency confound: high-pass cutoff (secs), Inf = no filtering, default is 128
nBase = 12;


% adjust the number of processors to work with to number of subjects------
if size(allSubs) < nProc
    nProc = size(allSubs,2);
end

% assign subjects to processors-------------------------------------------
% NOTE: splitSubs contain 1 x n_proc cells, each containing 1 or more subjects
splitSubs  = splitvect(allSubs, nProc);

%% define what batches to run
fprintf('\n##############################\nplease check your order now:\n##############################\n');

check = 0  % 0 to run, 1 to not run chosen batch

smoothie                = 0 %in case of smoothing realigned niftis: this has to be done first (separately); already done in previous step
disp('##############################');
specifyFirstLevel       = 1
estimateFirstLevel      = 1
estimateContrasts       = 1

if strcmp(modelname,'complete_onset_stimulus_condition_contrast')
    estimateContrasts_condition = 1
else
    estimateContrasts_condition = 0;
end

use_sNIFs               = 0 %in case of using the smoothed realigned niftis

disp('###############################');
% this has to be done separately
prepNormAndOrSmooth     = 0 %need this to do smooth and nl_smooth

%from UB-----------------------------------
norm2mni                = 0 %for now using those below from CB

%from CB-----------------------------------
nl_warp                 = 0
nl_smooth               = 0 % smooth only works when done separately

%%===================================================
%% loop across subjects
%%====================================================
for np = 1:size(splitSubs,2)

    matlabbatch = [];

    % loop across subjects whitin each splitSubs-cell---------------------
    for g = 1:size(splitSubs{np},2)

        % define subjectspecific values-----------------------------------
        subName = sprintf('sub-%02.2d',splitSubs{np}(g));
        epis = struct('name',subName);
        hdrDir  = [derivDir];
        fprintf(['\n\nDoing volunteer ' subName '\n']);

        % read in values from dicom header
        load(fullfile(hdrDir,['dcmhdr.mat']),'dcmhdr');
        tr              = dcmhdr.acqpar.RepetitionTime;
        trInSecs        = tr/1e3;
        pixelSpacing    = dcmhdr.acqpar.PixelSpacing;
        sliceThickness  = dcmhdr.acqpar.SliceThickness;
        vox             = [pixelSpacing' sliceThickness];
        VboundingBox    = [NaN NaN NaN; NaN NaN NaN];

        % temporal adjustments----------------------------------------------------
        adjValue   = 0;%1;%
        tempAdj    = adjValue*trInSecs;%move onsets (0 or adjustment value)
        baseRes    = nBase*trInSecs;

        %firstlevel name----------------------------------------------------------
        firstLevelName  = sprintf('%s_%s',basisF,modelname);
        % firstLevelName = 'test';
        fprintf('\nbasis function is %s\norto is %d\ntempAdj is %d\n',basisF,orthog,tempAdj);

        %create empty variables-------------------------------------------
        epiDir     = cell(size(allRuns,2)*2,1);
        epiFiles   = cell(size(epiDir,1),0);

        % Pharma condition Load in
        pharm_conds = '/home/nold/Desktop/PEEP/fMRI/Analysis/fmri_analysis/first_level/PEEP_entblindungsliste_dummy_coded_final.csv';
        [~,~, ~, treatment_order]= textread(pharm_conds, '%d %d %d %s', 'headerlines', 1, 'delimiter',',');
        treatment_order = str2num(cat(1,treatment_order{:}));

        %looping across sessions to get the files per session-----------------
        for ses = 1:length(session)

            %looping across runs to get the files per run---------------------
            for nr = 1:size(allRuns,2)
                runName        = sprintf('run-%02.2d',nr);
                strnameRun     = sprintf('ses_%02d_run_%02.2d',session(ses),nr);
                epiDir         = fullfile(derivDir, subName, subDir{2});
                anatDir        = fullfile(derivDir, subName, subDir{1});

                %% Specify first level
                %%=========================================================
                cd(epiDir);%go to location of epi-files

                if use_sNIFs == 1
                    epiName = ['s' num2str(skern) '_ra' subName sprintf('_ses-%02.2d',session(ses)) '_task-' studyName '_' runName '_bold.nii'];
                else
                    epiName = ['bra' subName sprintf('_ses-%02.2d',session(ses)) '_task-' studyName '_' runName '_bold.nii'];
                end

                numberScans = size(spm_select('expand',fullfile(epiDir,epiName)),1);


                %get image files
                %==========================================================
                if strcmp(imageDim, '4D')
                    epiFiles = cell(numberScans,1);

                    for j = 1:numberScans
                        epiFiles{j}   = [epiDir filesep epiName ',' num2str(j)];
                    end
                    clear j

                elseif strcmp(imageDim, '3D') %@UB: not yet tested this version


                end

                if isempty(epiFiles)
                    fprintf('Attention: NO epi files!');
                    return;
                end

                fprintf('\n%s has %d files in %s in session %d\n',subName,length(epiFiles),runName,session(ses));
                epis.(strnameRun) = epiFiles;%each epi-name-string is one cell



                %% Get movement and physio files (noise regressors)
                % Noise regressoers:
                % - 24 Movement regressors (derivatives)
                % - x Physio (pulse and respiration)
                % - x Spike (spike movements within sessions)
                %%=======================================================
                physio_dir      = '/projects/crunchie/nold/PEEP/physiology/MAIN/noise_regressors';
                regressors.(strnameRun) = fullfile(physio_dir, ['BS_complete_reg_' subName sprintf('_ses_%02.2d',session(ses)) '_' sprintf('run_%02d',allRuns(nr)) '.txt']);

            end %for the loop across runs

        end % end for loop across sessions

        clear nr runName;

        %% Setup epis for loop in matlabbatch
        %=======================================================
        if size(allRuns,2) == 1
            epis.run = {epis.run_01};
            regressors.runs = {regressors.run_01};
        else
            epis.run = [{epis.ses_01_run_01} {epis.ses_01_run_02} {epis.ses_01_run_03} {epis.ses_01_run_04} {epis.ses_02_run_01} {epis.ses_02_run_02} {epis.ses_02_run_03} {epis.ses_02_run_04}];
            regressors.runs = [{regressors.ses_01_run_01} {regressors.ses_01_run_02} {regressors.ses_01_run_03} {regressors.ses_01_run_04} {regressors.ses_02_run_01} {regressors.ses_02_run_02} {regressors.ses_02_run_03} {regressors.ses_02_run_04}];

        end



        %% set up 1level directory
        %==========================================================
        %define name of firstlevel-dir
        firstLevelDir = [FFXDir subName filesep];

        %set index for matlabbatch
        gi   = 1;

        if smoothie == 1

            % smooth niftis BEFORE specifying firstlevel
            %=============================================================
            %get realigned niftis-----------------------------------------
            for nr = 1:size(allRuns,2) %@UB: the loop works when doing one sub at a time. I might remove this section from the firstlevel script
                runName     = sprintf('run-%02d',nr);
                funcDir     = fullfile(derivDir, 'spm_preprocessing', subName, subDir{2});
                ntfDir      = funcDir; %[funcDir filesep runName filesep];
                rNIFs   = char(spm_select('ExtFPList', ntfDir, ['^bra' subName sprintf('_ses-%02.2d',session) '_task-' studyName '_' runName '_bold.nii'])); % slicetimed and realigned niftis (rasub-01_task-pdev_run-01_bold.nii)


                matlabbatch{gi}.spm.spatial.smooth.data   	= cellstr(rNIFs);
                matlabbatch{gi}.spm.spatial.smooth.fwhm   	= skernel;
                matlabbatch{gi}.spm.spatial.smooth.dtype  	= 0;
                matlabbatch{gi}.spm.spatial.smooth.im       = 0;
                matlabbatch{gi}.spm.spatial.smooth.prefix 	= ['s',num2str(skern) '_'];

                gi = gi + 1;
            end
        end %for smoothing

        %######################################################################################
        if specifyFirstLevel == 1

            fprintf('\nspecifying first level\nfirstLevelName is %s\n\n',firstLevelName);
            fprintf('\nCAUTION! If you continue, all content in folder %s will be deleted!!\n',firstLevelDir);
            fprintf('\nto continue press F5. To cancel press ctrl+C\n###################################\n');

            %=============================================================================
            if isfolder(firstLevelDir)
                rmdir((firstLevelDir),'s') %in case of old model ('s' includes subdir-tree)
            end

            mkdir(firstLevelDir)
            cd(firstLevelDir)
            clear SPM %in case of old model

            %firstlevel specify matlabbatch
            %=========================================================================
            matlabbatch{gi}.spm.stats.fmri_spec.dir = {firstLevelDir};
            %timing
            matlabbatch{gi}.spm.stats.fmri_spec.timing.units = units;
            matlabbatch{gi}.spm.stats.fmri_spec.timing.RT = trInSecs;%time of repetition (TR)
            matlabbatch{gi}.spm.stats.fmri_spec.timing.fmri_t = 16;%microtime resolution
            matlabbatch{gi}.spm.stats.fmri_spec.timing.fmri_t0 = 8;%microtime onset


            %% loop across runs for both session
            l = 1;
            for ses = 1:length(session)
                for nr = 1:size(allRuns,2)

                    runName = sprintf('run-%02d',nr);
                    temp = ['bra' subName sprintf('_ses-%02.2d',session(ses)) '_task-' studyName '_' runName '_bold.nii'];

                    matlabbatch{gi}.spm.stats.fmri_spec.sess(l).scans = cellstr(spm_select('ExtFPList',epiDir,temp));

                    %% read onsets----------------------------------------------
                    %NOTE:
                    %conds is a (1 x nConds) struct array with the fields 'name','onset', 'duration','pmod','orth'
                    %condsFIR is a (1 x nConds) struct array with the fields
                    %'name', 'onset', 'duration'; duration is defined as 0
                    % NOTE:

                    if strcmp(modelname,'complete_different_onsets_heat_pressure')
                        % DIFFERENT ONSETS Heat Onset Plateau Start, Pressure Onset Stimulus Start
                        [conds] = extract_onsets_peep_different_HEAT_and_PRESSURE(splitSubs{np}(g),session(ses),studyPart,nr,basisF);

                    elseif strcmp(modelname,'complete_onset_stimulus')
                        % SAME ONSETS Heat and ressure Onset: Stimulus Start
                        [conds] = extract_onsets_peep_ONSET_START_STIMULUS_HEAT_and_PRESSURE(splitSubs{np}(g),session(ses),nr,basisF);

                    elseif strcmp(modelname,'split_onset_stimulus')
                        % SPLIT HRF: EARLY AND LATE PAIN (ONSETS Heat and ressure Onset: Stimulus Start)
                        [conds] = extract_onsets_peep_ONSET_START_STIMULUS_SPLIT_ONSET(splitSubs{np}(g),session(ses),nr,basisF);

                    elseif strcmp(modelname,'complete_onset_stimulus_condition_contrast')
                        % SAME ONSETS Heat and ressure Onset: Stimulus Start
                        [conds] = extract_onsets_peep_ONSET_START_STIMULUS_HEAT_and_PRESSURE(splitSubs{np}(g),session(ses),nr,basisF);

                    end

                    matlabbatch{gi}.spm.stats.fmri_spec.sess(l).cond = conds;
                    matlabbatch{gi}.spm.stats.fmri_spec.sess(l).regress = struct('name', {}, 'val', {});
                    matlabbatch{gi}.spm.stats.fmri_spec.sess(l).multi_reg = {regressors.runs{l}};
                    matlabbatch{gi}.spm.stats.fmri_spec.sess(l).hpf = hpf;

                    l = l + 1;

                    clear conds condsFIR;

                end %end for loop nr across runs
            end
            clear l

            matlabbatch{gi}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});


            if strcmp(basisF, 'HRF')
                matlabbatch{gi}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];%appears in design as bf (basis function)
            elseif strcmp(basisF, 'FIR')
                matlabbatch{gi}.spm.stats.fmri_spec.bases.fir.length = baseRes; %nBase*TRInSec FIR model
                matlabbatch{gi}.spm.stats.fmri_spec.bases.fir.order = nBase; %FIR model
            end

            %NOTE:length is the length in seconds of the post-stimulus time window that the basis functions span.
            %order specifies how many basis functions to use.
            %In the FIR model we create one bin per maxStimLength in secs, each covering the time of one trInSecs.
            funcDir     = fullfile(derivDir, subName, subDir{2});
            coreg_skull_strip = sprintf('bins3cbrainstem_mask.nii');
            matlabbatch{gi}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{gi}.spm.stats.fmri_spec.global = 'None'; %global normalization: OPTINS:'Scaling'|'None'
            matlabbatch{gi}.spm.stats.fmri_spec.mthresh = -Inf;
            matlabbatch{gi}.spm.stats.fmri_spec.mask = {[funcDir filesep coreg_skull_strip]}; %smoothed and coregisteres skull-strip
            matlabbatch{gi}.spm.stats.fmri_spec.cvi = cvi;

            gi = gi + 1;
        end %end for specifyFirstLevel


        if estimateFirstLevel == 1
            %firstlevel estimate matlabbatch
            %=========================================================================
            matlabbatch{gi}.spm.stats.fmri_est.spmmat = {fullfile(firstLevelDir,'SPM.mat')};
            matlabbatch{gi}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{gi}.spm.stats.fmri_est.method.Classical = 1;

            gi = gi + 1;
        end %end for estimateFirstLevel

        if estimateContrasts == 1
            % Contrast Estimation
            %=========================================================================
            %NOTE: sessrep can take the values:
            %'none' do not replicate
            %'repl' replicate over sessions
            %'replsc' replicate and scale
            %'sess' create per session
            %'both' replicate and create per session
            %'bothsc' replicate and scale and create per session


            %% Extract number of noise regressors for later contrast definition
            % =========================================================
            clear n
            n_nuis_file = sprintf('n_nuis_sub-%02d',splitSubs{np}(g));
            n = load(fullfile(physio_dir,n_nuis_file));
            %n.n_nuis = zeros(1,8);

            template = [];
            template.spm.stats.con.spmmat = {fullfile(firstLevelDir,'SPM.mat')};
            template.spm.stats.con.delete = 1;

            %% Contrasts FIR

            if strcmp(basisF,'FIR')


                tco = 0;

                %T-con=====================================================
                %prep values-------------------------------------------
                nSess = 8;
                n_ints = nBase;
                co = 1;


                %% Heat 30 Nacl
                %====================

                for i_fir = 1:n_ints



                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Heat30_' num2str(i_fir) '_nacl'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end


                %% Heat 50 Nacl
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Heat50_' num2str(i_fir) '_nacl'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end



                %% Heat 70 Nacl
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Heat70_' num2str(i_fir) '_nacl'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*2 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end


                %% Pressure 30 Nacl
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Pressure30_' num2str(i_fir) '_nacl'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*3 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end



                %% Pressure 50 Nacl
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Pressure50_' num2str(i_fir) '_nacl'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*4 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end




                %% Pressure 70 Nacl
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Pressure70_' num2str(i_fir) '_nacl'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*5 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end





                %% Heat 30 Nax
                %====================

                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Heat30_' num2str(i_fir) '_nax'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);


                        if treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end


                %% Heat 50 Nax
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Heat50_' num2str(i_fir) '_nax'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints + i_fir) = 1;

                    for j = 1:nSess


                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end



                %% Heat 70 Nax
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Heat70_' num2str(i_fir) '_nax'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*2 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);
                        if treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end

                %% Pressure 30 Nax
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Pressure30_' num2str(i_fir) '_nax'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*3 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end



                %% Pressure 50 Nax
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Pressure50_' num2str(i_fir) '_nax'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*4 + i_fir) = 1;

                    for j = 1:nSess

                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end




                %% Pressure 70 Nax
                %=========================
                for i_fir = 1:n_ints

                    tco = tco + 1;
                    con = [];
                    template.spm.stats.con.consess{tco}.tcon.name   = ['Pressure70_' num2str(i_fir) '_nax'];
                    tpl1        = zeros(1,n_ints*6);
                    tpl2        = zeros(1,n_ints*6);
                    tpl1(n_ints*5 + i_fir) = 1;

                    for j = 1:nSess
                        nMp = n.n_nuis(j);
                        mvmnt = zeros(1,nMp);

                        if treatment_order(splitSubs{np}(g)) == 1

                            if j < floor(nSess/2)+1
                                con = [con tpl1 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl2 mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0

                            if j < floor(nSess/2)+1
                                con = [con tpl2 mvmnt];
                            elseif j > floor(nSess/2)
                                con = [con tpl1 mvmnt];
                            end

                        end

                    end
                    template.spm.stats.con.consess{tco}.tcon.convec  = con;
                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    co = co + 1; %increment by 1
                end



               
                %% Contrasts HRF

            elseif strcmp(basisF, 'HRF') && strcmp(modelname,'complete_onset_stimulus')
                if ~estimateContrasts_condition
                    tco = 0;
                    % ===================================================
                    %T-cons
                    % ==================================================

                    % General Values

                    nConds = 6;

                    %Pain Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pain_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length          
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                         error('Contrast Length is not coherent');
                    end


                    % Pain Nax
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pain_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Heat Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                
                    % Heat Nax
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end
                    


                    %Pressure Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Pressure Nax
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))+length(subjectConstant)
                        error('Contrast Length is not coherent');
                    end


                                        

                     
                    % Heat NAX > Nacl
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                    subjectConstant = [];%zeros(1,nSess);
                    zeroBlock = [zeros(1,nConds) 0];

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nlx>sal'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    % Heat SAL > NLX
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                    subjectConstant = [];%zeros(1,nSess);
                    zeroBlock = [zeros(1,nConds) 0];

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_sal>nlx'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    % Pressure NAX > Nacl
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    negcondBlock = [zeros(1,nConds/2)  -ones(1,nConds/2) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nlx>sal'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Pressure NAX > Nacl
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    negcondBlock = [zeros(1,nConds/2)  -ones(1,nConds/2) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_sal>nlx'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negcondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                                        
                    % Heat SAL  > Pressure SAL
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) -ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_sal>pressure_sal'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end




                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Pressure SAL  > Heat SAL
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [-ones(1,nConds/2) ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_sal>heat_sal'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end




                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                   
                    % Heat NLX  > Pressure NLX
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) -ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nlx>pressure_nlx'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end




                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Pressure NLX  > Heat NLX
                    %==================
                                         
                    % Prep values
                    tConVec = [];
                    condBlock = [-ones(1,nConds/2) ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nlx>heat_nlx'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end




                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                   %% Parametric effect
                                      
                    %% Heat SAL
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_int_sal'};

                         %convec-------------------------------------------------
                         condBlock = [-1 0 1 zeros(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    %% Heat NLX
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);
                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_int_nlx'};

                         %convec-------------------------------------------------
                         condBlock = [-1 0 1 zeros(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    %% Pressure SAL
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_int_sal'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,nConds/2) -1 0 1 0];
                         zeroBlock = [zeros(1,nConds) 0];

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    %% Pressure NLX
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);
                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_int_nlx'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,nConds/2) -1 0 1 0];
                         zeroBlock = [zeros(1,nConds) 0];

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     %% Exercise

                    %% Heat Exercise Nacl
                    %==================
                      
                    % Ex low

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo_sal'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock;
                        elseif exercise(r) == 1
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Ex High

                      % Ex low

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exHi_sal'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 1
                            condBlock;
                        elseif exercise(r) == 0
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                                        
                    %% Heat Exercise NLX
                    %==================
                      
                    % Ex low

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo_nlx'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock;
                        elseif exercise(r) == 1
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Ex High

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exHi_nlx'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 1
                            condBlock;
                        elseif exercise(r) == 0
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Pressure Exercise Nacl
                    %==================
                      
                    % Ex low

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exLo_sal'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock;
                        elseif exercise(r) == 1
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Ex High
                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exHi_sal'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 1
                            condBlock;
                        elseif exercise(r) == 0
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                                        
                    %% Pressure Exercise NLX
                    %==================
                      
                    % Ex low

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exLo_nlx'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock;
                        elseif exercise(r) == 1
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Ex High

                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exHi_nlx'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 1
                            condBlock;
                        elseif exercise(r) == 0
                            condBlock = zeroBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% Exercise INduced Hypolagesia (Low vs. High Intensity)
                    
                    %% Heat Nacl

                    % Heat Exercise High > Low Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exHi>exLo_sal'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = -condBlock;
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 1
                            condBlock;
                        elseif exercise(r) == 0
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Heat Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo>exHi_sal'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = -condBlock;
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% Pressure Nacl                   
                    
                    % Pressure Exercise High > Low Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exHi>exLo_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                         negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    % Pressure Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exLo>exHi_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                         negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock= negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    %%  Heat NLX

                    % Heat Exercise High > Low Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exHi>exLo_nlx'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Heat Exercise Low > High NLX
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
                    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo>exHi_nlx'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [ones(1,nConds/2) zeros(1,nConds/2) 0];
                        negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% Pressure NLX                   
                    
                    % Pressure Exercise High > Low NLX
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exHi>exLo_nlx'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                         negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    % Pressure Exercise Low > High NLX
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                    negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_exLo>exHi_nlx'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
                         negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock= negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    

                    %% EIH and Intensities

                    %% Heat 70 SAL
                    % Heat Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat70_exHi>exLo_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [0 0 1 zeros(1,nConds/2) 0];
                         negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    condBlock = [0 0 1 zeros(1,nConds/2) 0];
                    negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat70_exLo>exHi_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [0 0 1 zeros(1,nConds/2) 0];
                         negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Heat 70 NLX
                    % Heat Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat70_exHi>exLo_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [0 0 1 zeros(1,nConds/2) 0];
                        negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat70_exLo>exHi_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [0 0 1 zeros(1,nConds/2) 0];
                        negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end



                    %% Heat 50 SAL
                    % Heat Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat50_exHi>exLo_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [0 1 0 zeros(1,nConds/2) 0];
                         negcondBlock = [0 -1 0 zeros(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];

                         if exercise(r) == 0
                             condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat50_exLo>exHi_sal'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [0 1 0 zeros(1,nConds/2) 0];
                        negcondBlock = [0 -1 0 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock;
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Heat 50 NLX
                    % Heat Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat50_exHi>exLo_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [0 1 0 zeros(1,nConds/2) 0];
                        negcondBlock = [0 -1 0 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat50_exLo>exHi_nlx'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [0 1 0 zeros(1,nConds/2) 0];
                         negcondBlock = [0 -1 0 zeros(1,nConds/2) 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                                        
                    %% Heat 30 SAL
                    % Heat Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat30_exHi>exLo_sal'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [1 0 0 zeros(1,nConds/2) 0];
                        negcondBlock = [-1 0 0 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat30_exLo>exHi_sal'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [1 0 0 zeros(1,nConds/2) 0];
                        negcondBlock = [-1 0 0 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Heat 30 NLX
                    % Heat Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat30_exHi>exLo_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [1 0 0 zeros(1,nConds/2) 0];
                        negcondBlock = [-1 0 0 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                 
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat30_exLo>exHi_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [1 0 0 zeros(1,nConds/2) 0];
                        negcondBlock = [-1 0 0 zeros(1,nConds/2) 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                       
                    %% Pressure 70 SAL
                    % Pressure Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure70_exHi>exLo_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) 0 0 1 0];
                         negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure70_exLo>exHi_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) 0 0 1 0];
                         negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
                         zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Pressure 70 NLX
                    % Pressure Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure70_exHi>exLo_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 0 0 1 0];
                        negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure70_exLo>exHi_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 0 0 1 0];
                        negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end



                    %% Pressure 50 SAL
                    % Pressure Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure50_exHi>exLo_sal'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) 0 1 0 0];
                         negcondBlock = [zeros(1,nConds/2) 0 -1 0 0];
                         zeroBlock = [zeros(1,nConds) 0];

                         if exercise(r) == 0
                             condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure50_exLo>exHi_sal'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 0 1 0 0];
                        negcondBlock = [zeros(1,nConds/2) 0 -1 0 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock;
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Pressure 50 NLX
                    % Pressure Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure50_exHi>exLo_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 0 1 0 0];
                        negcondBlock = [zeros(1,nConds/2) 0 -1 0 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure50_exLo>exHi_nlx'};

                         %convec-------------------------------------------------
                         exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                         condBlock = [zeros(1,nConds/2) 0 1 0 0];
                         negcondBlock = [zeros(1,nConds/2) 0 -1 0 0];
                         zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                                        
                    %% Pressure 30 SAL
                    % pressure Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure30_exHi>exLo_sal'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 1 0 0 0];
                        negcondBlock = [zeros(1,nConds/2) -1 0 0 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High Nacl
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                   
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure30_exLo>exHi_sal'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 1 0 0 0];
                        negcondBlock = [zeros(1,nConds/2) -1 0 0 0];
                        zeroBlock = [zeros(1,nConds) 0];

                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Pressure 30 NLX
                    % Pressure Exercise High > Low 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure30_exHi>exLo_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 1 0 0 0];
                        negcondBlock = [zeros(1,nConds/2) -1 0 0 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock = negcondBlock;
                        elseif exercise(r) == 1
                            condBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Exercise Low > High 
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                 
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure30_exLo>exHi_nlx'};

                        %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock = [zeros(1,nConds/2) 1 0 0 0];
                        negcondBlock = [zeros(1,nConds/2) -1 0 0 0];
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    % Interaction MODALITY x TREATMENT
                    % (Heat NLX > Heat SAL) > (Pressure NLX > Pressure SAL)
                    % ===================

                                           
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/2) -ones(1,nConds/2) 0];
                    negCondBlock = [-ones(1,nConds/2) ones(1,nConds/2) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_modality_treatment'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end



                    % Interaction INTENSITY x TREATMENT
                    % Int NLX > INT SAL
                    % ===================

                                           
                    % Prep values
                    tConVec = [];
                    condBlock = [-1 0 1 -1 0 1 0];
                    negCondBlock = [1 0 -1 1 0 -1 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Interaction INTENSITY x TREATMENT HEAT
                    % ===================
               
                    % Prep values
                    tConVec = [];
                    condBlock = [-1 0 1 0 0 0 0];
                    negCondBlock = [1 0 -1 0 0 0 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment_heat'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Interaction INTENSITY x TREATMENT Pressure
                    % ===================
               
                    % Prep values
                    tConVec = [];
                    condBlock = [0 0 0 -1 0 1 0];
                    negCondBlock = [0 0 0 1 0 -1 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment_pressure'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Interaction INTENSITY x MODALITY
                    % ===================

                                           
                    % Prep values
                    tConVec = [];
                    condBlock = [-1 0 1 1 0 -1 0];
                    %negCondBlock = [1 0 -1 1 0 -1 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_modality'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

%                         
%                     % Interaction Exercise x MODALITY (HEAT)
%                     % (Heat ExLoNLX > Heat ExHiNLX) > (Heat ExLoSAL > Heat ExHiSAL)
%                     % ===================
% 
%                                            
%                     % Prep values
%                     tConVec = [];
%                     condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
%                     negCondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
%                     subjectConstant = [];%zeros(1,nSess);
% 
%                     
%                     for r = 1:nSess
% 
%                         nMp = n.n_nuis(r);
%                         mvmnt = zeros(1,nMp);
% 
%                         %name--------------------------------------------------
%                         tcoName = {'inter_exercise_treatment'};
% 
%                          %convec------------------------------------------------
% 
%                          if r < floor(nSess/2)+1
%                              if exercise(r) == 0
%                                  condBlock;
%                              elseif exercise(r) == 1
%                                  condBlock = negCondBlock;
%                              end
%                          elseif r < floor(nSess/2)+1
%                               if exercise(r) == 0
%                                  condBlock = negCondBlock;
%                              elseif exercise(r) == 1
%                                  condBlock;
%                              end
%                          end
% 
% 
% 
%                         if treatment_order(splitSubs{np}(g)) == 1
%                             if r < floor(nSess/2)+1
%                                 tConVec = [tConVec condBlock mvmnt subjectConstant];
%                             elseif r > floor(nSess/2)
%                                 tConVec = [tConVec negCondBlock mvmnt subjectConstant];
%                             end
% 
%                         elseif treatment_order(splitSubs{np}(g)) == 0
%                             if r < floor(nSess/2)+1
%                                 tConVec = [tConVec negCondBlock mvmnt subjectConstant];
%                             elseif r > floor(nSess/2)
%                                 tConVec = [tConVec condBlock mvmnt subjectConstant];
%                             end
%                         end
% 
%                     end
% 
%                     %fill t-cons into template
%                     for co = 1:size(tcoName,2)
%                         tco = tco + 1;
%                         template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
%                         template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
%                         template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
%                     end
% 
%                     % Sanity Check: Contrast Vector length
%                     if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
%                         error('Contrast Length is not coherent');
%                     end


  
                    % Heat 30 Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [1 zeros(1,nConds-1) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat30_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Heat 50 Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [0 1 0 zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat50_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Heat 70 Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [0 0 1 zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat70_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Pressure 30 Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) 1 0 0 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure30_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Pressure 50 Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) 0 1 0 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure50_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Pressure 70 Nacl
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) 0 0 1 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure70_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                      % Heat 30 nlx
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [1 zeros(1,nConds-1) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat30_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Heat 50 NLX
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [0 1 0 zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat50_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Heat 70 NLX
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [0 0 1 zeros(1,nConds/2) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat70_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                        
                    % Pressure 30 NLX
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) 1 0 0 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure30_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Pressure 50 NLX
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) 0 1 0 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure50_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     % Pressure 70 NLX
                    %==================

                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds/2) 0 0 1 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure70_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    % Rating
                    %==================
            
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,nConds) 1];
                    zeroBlock = [zeros(1,nConds) 1];
                    subjectConstant = [];%zeros(1,nSess);
                   

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'rating'};

                        %convec-------------------------------------------------
                        tConVec = [tConVec condBlock mvmnt subjectConstant];

                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end)) + length(subjectConstant)
                        error('Contrast Length is not coherent');
                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end
                 end

            elseif strcmp(basisF, 'HRF') && strcmp(modelname,'split_onset_stimulus')

                    nConds = 12;
                    nSess = 8;
                    tco = 0;

                    %Heat Nacl Early Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/4) zeros(1,nConds-nConds/4) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nacl_early'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                     %Heat Nacl Late Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) ones(1,3) zeros(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nacl_late'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                                        
                    %Pressire Nacl Early Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,3),ones(1,3) zeros(1,6) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nacl_early'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                     %Heat Nacl Late Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,9) ones(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nacl_late'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                     %Heat nlx Early Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [ones(1,nConds/4) zeros(1,nConds-nConds/4) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nlx_early'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Heat nlx Late Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) ones(1,3) zeros(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_nlx_late'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                                        
                    %Pressure nlx Early Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,3),ones(1,3) zeros(1,6) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nlx_early'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    %Pressure nlx Late Pain
                    %=================                    
                
                    % Pressure nlx Late Pain
                    %=================                    
                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,9) ones(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_nlx_late'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    %% Early Pain Heat

                    %% NLX > SAL

                    %Prep values
                    tConVec = [];
                    condBlock = [ones(1,3) zeros(1,9) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_early_nlx_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end



                    %% SAL > NLX

                     %Prep values
                    tConVec = [];
                    condBlock = [ones(1,3) zeros(1,9) 0];
                    negcondBlock = [-ones(1,3) zeros(1,9) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_early_sal_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end




                    %% Late Pain Heat

                    %% NLX > SAL

                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) ones(1,3) zeros(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_late_nlx_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    %% SAL > NLX

                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) ones(1,3) zeros(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_late_sal_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end



                    %% Early Pain Pressure

                    %% NLX > SAL

                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,3),ones(1,3) zeros(1,6) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_early_nlx_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end



                    %% SAL > NLX

                     %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,3),ones(1,3) zeros(1,6) 0];
                    negcondBlock = [-ones(1,3) zeros(1,9) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_early_sal_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end




                    %% Late Pain Pressure

                    %% NLX > SAL

                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,9) ones(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_late_nlx_sal'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end
 
 
                     %% SAL > NLX

                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,9) ones(1,3) 0];
                    zeroBlock = [zeros(1,nConds) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_late_sal_nlx'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end



                    %% Parametric effect
                                      
                    %% Heat INT SAL EARLY
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_int_sal_early'};

                         %convec-------------------------------------------------
                         condBlock = [-1 0 1 zeros(1,9) 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                  

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                      
                    %% Heat INT SAL LATE
                    %==================
                                                            
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_int_sal_late'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,6) -1 0 1 zeros(1,3) 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                  

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    %% Heat INT NLX EARLY
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);
                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_int_nlx_early'};

                         %convec-------------------------------------------------
                         condBlock = [-1 0 1 zeros(1,9) 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end



                    %% Heat INT NLX LATE
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);
                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_int_nlx_late'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,6) -1 0 1 zeros(1,3) 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                                        
                    %% PRESSURE INT SAL EARLY
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_int_sal_early'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,3) -1 0 1 zeros(1,6) 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                  

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                                      
                    %% PRESSURE INT SAL LATE
                    %==================
                                                            
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_int_sal_late'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,9) -1 0 1 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                  

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    %% PRESSURE INT NLX EARLY
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);
                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_int_nlx_early'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,3) -1 0 1 zeros(1,6) 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end



                    %% PRESSURE INT NLX LATE
                    %==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);
                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'pressure_int_nlx_late'};

                         %convec-------------------------------------------------
                         condBlock = [zeros(1,9) -1 0 1 0];
                         zeroBlock = [zeros(1,12) 0];

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% INTERACTION TREATMENT AND EARLY/LATE PAIN HEAT
                    %% ========================

                    %Prep values
                    tConVec = [];
                    condBlock = [-ones(1,3) zeros(1,3) ones(1,3) zeros(1,3) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_treatment_early_late_heat'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    %% INTERACTION TREATMENT AND EARLY/LATE PAIN Pressure
                    % =============================

                    %Prep values
                    tConVec = [];
                    condBlock = [zeros(1,3)-ones(1,3) zeros(1,3) ones(1,3) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_treatment_early_late_pressure'};

                        %convec-------------------------------------------------
                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec -condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Interaction INTENSITY x TREATMENT HEAT EARLY
                    % ===================
               
                    % Prep values
                    tConVec = [];
                    condBlock = [-1 0 1 0 0 0 zeros(1,6) 0];
                    negCondBlock = [1 0 -1 0 0 0 zeros(1,6) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment_heat_early'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                    % Interaction INTENSITY x TREATMENT HEAT LATE
                    % ===================
               
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) -1 0 1 0 0 0 0];
                    negCondBlock = [zeros(1,6) 1 0 -1 0 0 0 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment_heat_late'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                     % Interaction INTENSITY x TREATMENT Pressure EARLY
                    % ===================
               
                    % Prep values
                    tConVec = [];
                    condBlock = [0 0 0 -1 0 1 zeros(1,6) 0];
                    negCondBlock = [0 0 0 1 0 -1 zeros(1,6) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment_pressure_early'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                                    
                    % Interaction INTENSITY x TREATMENT Pressure LATE
                    % ===================
               
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) 0 0 0 -1 0 1 0];
                    negCondBlock = [zeros(1,6) 0 0 0 1 0 -1 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_intensity_treatment_pressure_late'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% Interaction MODALITY x TREATMENT EARLY
                    %% ===================

                                           
                    % Prep values
                    tConVec = [];
                    condBlock = [ones(1,3) -ones(1,3) zeros(1,6) 0];
                    negCondBlock = [-ones(1,3) ones(1,3) zeros(1,6) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_modality_treatment_early'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% Interaction MODALITY x TREATMENT LATE
                    %% ===================

                                           
                    % Prep values
                    tConVec = [];
                    condBlock = [zeros(1,6) ones(1,3) -ones(1,3) 0];
                    negCondBlock = [zeros(1,6) -ones(1,3) ones(1,3) 0];
                    subjectConstant = [];%zeros(1,nSess);

                    
                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'inter_modality_treatment_late'};

                         %convec------------------------------------------------

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec negCondBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end

                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                    %% Heat Exercise Low > High Nacl Early
                    %% ==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo>exHi_sal_early'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [ones(1,3) zeros(1,9) 0];
                        negcondBlock = -condBlock;
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     %% Heat Exercise Low > High Nacl Late
                    %% ==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo>exHi_sal_late'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [zeros(1,3) ones(1,3) zeros(1,6) 0];
                        negcondBlock = -condBlock;
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end




                    %% Heat Exercise Low > High NLX Early
                    %% ==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo>exHi_nlx_early'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [ones(1,3) zeros(1,9) 0];
                        negcondBlock = -condBlock;
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end

                     %% Heat Exercise Low > High NLX Late
                    %% ==================
                                        
                    % Prep values
                    tConVec = [];
                    subjectConstant = [];%zeros(1,nSess);

                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

                    for r = 1:nSess

                        nMp = n.n_nuis(r);
                        mvmnt = zeros(1,nMp);

                        %name--------------------------------------------------
                        tcoName = {'heat_exLo>exHi_nlx_late'};

                         %convec-------------------------------------------------
                        exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
                        condBlock =  [zeros(1,3) ones(1,3) zeros(1,6) 0];
                        negcondBlock = -condBlock;
                        zeroBlock = [zeros(1,nConds) 0];
                       
                        if exercise(r) == 0
                            condBlock; 
                        elseif exercise(r) == 1
                            condBlock = negcondBlock;
                        end

                        if treatment_order(splitSubs{np}(g)) == 1
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            end

                        elseif treatment_order(splitSubs{np}(g)) == 0
                            if r < floor(nSess/2)+1
                                tConVec = [tConVec zeroBlock mvmnt subjectConstant];
                            elseif r > floor(nSess/2)
                                tConVec = [tConVec condBlock mvmnt subjectConstant];
                            end
                        end

                    end

                    %fill t-cons into template
                    for co = 1:size(tcoName,2)
                        tco = tco + 1;
                        template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
                        template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                        template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
                    end


                    % Sanity Check: Contrast Vector length
                    if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
                        error('Contrast Length is not coherent');
                    end


                elseif estimateContrasts_condition && strcmp(modelname,'complete_onset_stimulus_condition_contrast')

                    %% Condition Contrast
                    mod = {'heat','pressure'};
                    int = {30,50,70};
                    ex = {'lo','hi'};
                    pharm = {'nacl','nax'};
                    nConds = 6;
                    extraCond = 1;
                    nRuns = 8; 
                    params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');
                    
                    heat_30         = [1 0 0 0 0 0 0];
                    heat_50         = [0 1 0 0 0 0 0];
                    heat_70         = [0 0 1 0 0 0 0];
                    pressure_30     = [0 0 0 1 0 0 0];
                    pressure_50     = [0 0 0 0 1 0 0];
                    pressure_70     = [0 0 0 0 0 1 0];
                    zeroBlock       = zeros(1,nConds+extraCond);
                    sessionConstant = [];
                    counter         = 1;
                    counterEx       = 1;
                    tco = 0;


                    for p = 1:numel(pharm)
                        for m = 1:numel(mod)
                            for e = 1:numel(ex)
                                for i = 1:numel(int)

                                    clear tConVec
                                    string = [mod{m} '_' num2str(int{i}) '_' ex{e} '_' pharm{p}];

                                    if strcmp(mod{m},'heat')
                                        if int{i} == 30
                                            condBlock = heat_30;
                                        elseif int{i} == 50
                                            condBlock = heat_50;
                                        elseif int{i} == 70
                                            condBlock = heat_70;
                                        end
                                    elseif strcmp(mod{m},'pressure')
                                        if int{i} == 30
                                            condBlock = pressure_30;
                                        elseif int{i} == 50
                                            condBlock = pressure_50;
                                        elseif int{i} == 70
                                            condBlock = pressure_70;
                                        end
                                    end

                                    %movement
                                    b1 = zeros(1,n.n_nuis(1));
                                    b2 = zeros(1,n.n_nuis(2));
                                    b3 = zeros(1,n.n_nuis(3));
                                    b4 = zeros(1,n.n_nuis(4));
                                    b5 = zeros(1,n.n_nuis(5));
                                    b6 = zeros(1,n.n_nuis(6));
                                    b7 = zeros(1,n.n_nuis(7));
                                    b8 = zeros(1,n.n_nuis(8));


                                    if counterEx < 4
                                        condBlock1 = condBlock;
                                        condBlock2 = zeroBlock;
                                    elseif counterEx > 3
                                        condBlock1 = zeroBlock;
                                        condBlock2 = condBlock;
                                    end

                                    counterEx = counterEx + 1;

                                    if counterEx > nConds
                                        counterEx = 1;
                                    end


                                    % extract exercise intensities
                                    if isequal(params_ses.P.exercise.condition,[0 0 1 1])
                                        conVec_ses1 = [condBlock1 b1 condBlock1 b2 condBlock2 b3 condBlock2 b4];
                                        conVec_ses2 = [condBlock1 b5 condBlock1 b6 condBlock2 b7 condBlock2 b8];
                                    elseif isequal(params_ses.P.exercise.condition,[0 1 1 0])
                                        conVec_ses1 = [condBlock1 b1 condBlock2 b2 condBlock2 b3 condBlock1 b4];
                                        conVec_ses2 = [condBlock1 b5 condBlock2 b6 condBlock2 b7 condBlock1 b8];
                                    elseif isequal(params_ses.P.exercise.condition,[1 1 0 0])
                                        conVec_ses1 = [condBlock2 b1 condBlock2 b2 condBlock1 b3 condBlock1 b4];
                                        conVec_ses2 = [condBlock2 b5 condBlock2 b6 condBlock1 b7 condBlock1 b8];
                                    elseif isequal(params_ses.P.exercise.condition,[0 1 0 1])
                                        conVec_ses1 = [condBlock1 b1 condBlock2 b2 condBlock1 b3 condBlock2 b4];
                                        conVec_ses2 = [condBlock1 b5 condBlock2 b6 condBlock1 b7 condBlock2 b8];
                                    elseif isequal(params_ses.P.exercise.condition,[1 0 1 0])
                                        conVec_ses1 = [condBlock2 b1 condBlock1 b2 condBlock2 b3 condBlock1 b4];
                                        conVec_ses2 = [condBlock2 b5 condBlock1 b6 condBlock2 b7 condBlock1 b8];
                                    elseif isequal(params_ses.P.exercise.condition,[1 0 0 1])
                                        conVec_ses1 = [condBlock2 b1 condBlock1 b2 condBlock1 b3 condBlock2 b4];
                                        conVec_ses2 = [condBlock2 b5 condBlock1 b6 condBlock1 b7 condBlock2 b8];
                                    end

                    
                                   
                                    if treatment_order(splitSubs{np}(g)) == 0
                                        if counter < 13
                                            sesZero = zeros(1,length(conVec_ses2));
                                            conVec = conVec_ses1;
                                            tConVec = [conVec sesZero sessionConstant];
                                        elseif counter > 12
                                            sesZero = zeros(1,length(conVec_ses1));
                                            conVec = conVec_ses2;
                                            tConVec = [sesZero conVec sessionConstant];
                                        end
                                    elseif treatment_order(splitSubs{np}(g)) == 1
                                        if counter < 13
                                            sesZero = zeros(1,length(conVec_ses1));
                                            conVec = conVec_ses2;
                                            tConVec = [sesZero conVec sessionConstant];
                                        elseif counter > 12
                                            sesZero = zeros(1,length(conVec_ses2));
                                            conVec = conVec_ses1;
                                            tConVec = [conVec sesZero sessionConstant];
                                        end
                                    end
                                
                                    counter = counter + 1;

                                    %fill into template
                                    tco = tco + 1;
                                    template.spm.stats.con.consess{tco}.tcon.name    = string;
                                    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
                                    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';


                                    % Sanity Check: Contrast Vector length
                                    if length(tConVec) ~= (length(conVec)+sessionConstant)*8 + length(sesZero)*8 +sum(n.n_nuis(1:end))
                                        error('Contrast Length is not coherent');
                                    end

                                end
                            end
                        end
                    end
                



            end %across basisF's (FIR or HRF)
   
            matlabbatch{gi} = template; %now add contrasts

            gi = gi + 1;
        end %for contrastEstimation

     




        %% run matlab (splitSubs)-----------------------------------
        if ~isempty(matlabbatch)
            fprintf('\nsplitSubs{np}(g) is %d\n',splitSubs{np}(g))
            spm_jobman('run',matlabbatch);
            %run_matlab(splitSubs{np}(g), matlabbatch, check);
        end

    end %for loop g across splitSubs


end %end for loop np across allSubs


%% function collection#####################################################

    function chuckCell = splitvect(v, n)
        % Splits a vector into number of n chunks of  the same size (if possible).
        % In not possible the chunks are almost of equal size.
        %
        % based on http://code.activestate.com/recipes/425044/
        chuckCell  = {};
        vectLength = numel(v);
        splitsize  = 1/n*vectLength;

        for i = 1:n
            idxs = [floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1]+1;
            chuckCell{end + 1} = v(idxs);
        end
    end

    function [files] = select_scans(file_filt, direc)
        x     = []; x     = spm_select('List', direc, file_filt);
        y     = []; y     = [repmat([direc filesep], size(x, 1), 1) x repmat(',1', size(x, 1), 1)];
        files = []; files = mat2cell(y, ones(size(y, 1), 1), size(y, 2));
    end
    function f_files = create_func_files(s_dir,f_templ,n_files)
        for i=1:n_files
            f_files{i} = [s_dir filesep f_templ ',' num2str(i)];
        end
    end
    function out = ins_letter(pscan,letter_start,letter_end)
        %inserts letter into name
        if nargin <3
            letter_end = [];
        end
        for a=1:size(pscan,1)
            [path , filename, extension] = fileparts(pscan(a,:));
            out(a,:) = [path filesep letter_start filename letter_end extension];
        end

    end

% function run_matlab(np, matlabbatch, check)
    function run_matlab(nsub, matlabbatch, check)

        spm_path          = fileparts(which('spm')); %get spm path
        mat_name          = which(mfilename);
        [~,mat_name,~]    = fileparts(mat_name);

        % fname = [mat_name '_' num2str(np) '.mat'];
        fname = [mat_name '_' num2str(nsub) '.mat'];

        save(fname,'matlabbatch');

        [~, name_stem]      = fileparts(fname);
        name_stem_m         = [name_stem '.m'];
        log_name            = strcat(name_stem, '.log');

        lo_cmd  = sprintf('clear matlabbatch;\nload(''%s'');\n',fname);
        ex_cmd  = sprintf('addpath(''%s'');\nspm(''defaults'',''FMRI'');\nspm_jobman(''initcfg'');\nspm_jobman(''run'',matlabbatch);\n',spm_path);
        end_cmd = sprintf('delete(''%s'');\ndelete(''%s'');\n',fname,name_stem_m);
        str     = strvcat(sprintf('function %s', name_stem),lo_cmd, ex_cmd, end_cmd, 'exit');

        spm_save(name_stem, str); %spm_save does not do .m files ...
        movefile(name_stem, name_stem_m); %rename

        if isunix
            matlab_exe  = [matlabroot filesep 'bin' filesep 'matlab']; % should be OK for unix
            matlab_lic  = license;
            cmd         = sprintf('xterm -e ''%s -nodesktop -nosplash -logfile %s -r "%s" '' &',matlab_exe, log_name, name_stem);
            %pointing to other license
            %     cmd         = sprintf('xterm -e ''%s -c "%s" -nodesktop -nosplash -logfile %s -r "%s" '' &',matlab_exe, matlab_lic, log_name, name_stem);
        elseif ispc
            matlab_exe  = 'matlab';
            cmd         = sprintf('start %s -nodesktop -nosplash -logfile %s -r "%s" exit',matlab_exe, log_name, name_stem);
        end

        if ~check
            system(cmd);
        end
    end
    function ons = extract_onsets(onset_file,logName)

        tf = load(onset_file);
        nevent = size(tf.p.log.events,1);
        nconds = size(logName,2);
        cnt = ones(1,nconds);

        for k = 1:nevent
            entry = tf.p.log.events(k,5);%number event, column 5
            for m = 1:nconds
                if iscell(entry{1}(2))
                    if strcmp(entry{1}(2),logName(m))
                        ons(m).onset(cnt(1,m)) = tf.p.log.events{k,3};
                        cnt(1,m) = cnt(1,m) + 1;
                    end
                elseif strcmp(entry{1},logName(m))
                    ons(m).onset(cnt(1,m)) = tf.p.log.events{k,3};
                    cnt(1,m) = cnt(1,m) + 1;
                else; continue;
                end
            end
        end
    end
%the original extract from CB----------------------------------------------
%reads out onsets in seconds
    function RES = extract_offset(onset_file,conditions)
        i1=1;i2=1;i3=1;i4=1;
        r1 = load(onset_file);
        nl = size(r1.p.log.events,1);
        for j=1:nl
            entry = r1.p.log.events(j,5);
            if iscell(entry{1}(2))
                switch entry{1}{2}
                    case conditions{1}
                        RES{1}.onset(i1) = r1.p.log.events{j,3};i1=i1+1;
                    case conditions{2}
                        RES{2}.onset(i2) = r1.p.log.events{j,3};i2=i2+1;
                    case conditions{3}
                        RES{3}.onset(i3) = r1.p.log.events{j,3};i3=i3+1;
                    case conditions{4}
                        RES{4}.onset(i4) = r1.p.log.events{j,3};i4=i4+1;
                end
            end
        end
        for i=1:size(conditions,2)
            RES{i}.name  = conditions{i};
        end
        % finally take 1st and last ss ff as rating trials
        RES{5}.onset = RES{1}.onset([1 end]);RES{1}.onset([1 end]) = [];
        RES{6}.onset = RES{2}.onset([1 end]);RES{2}.onset([1 end]) = [];
    end
end
%% ToDo-------------------------------------------------------------------
% warp and smooth do not work in one go
