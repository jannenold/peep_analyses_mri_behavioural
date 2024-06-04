function smooth_peep

addpath(genpath('/home/nold/Desktop/PEEP/fMRI/'));

%% call hostname and define paths and number of processing units
hostName    = char(getHostName(java.net.InetAddress.getLocalHost));
studyName   = 'peep';
session = 1;
studyPart = 'Pilot-01';
switch hostName
    case 'revelations.nin.uke.uni-hamburg.de'
        baseDir          = ['/projects/crunchie/nold/PEEP/fMRI/Data',filesep,studyPart,filesep,'rawdata/']; % übergeordnet
        derivDir         = ['/projects/crunchie/nold/PEEP/fMRI/Data',filesep,studyPart,filesep, 'derivatives/']; % untergeordnet
        nProc            = 8; % 8 parallel processes on revelations
        subDir            = {'anat','func','fmap'};
      
        
       
    otherwise
        error('host not recognized');
end

%% define general studyspecific values / names

%subjects / runs / conditions
%NOTE: excluded subs [1 6 25]; included subs [2:5 7:24 26:33];
allSubs     = 2;%[2:5 7:14];%26:33;% %it always processes max. 12
allRuns     = 1:4;
nSess       = numel(allRuns);
nCond       = 2;
condName    = {'pressure_on'; 'heat_on'};
conditions  = {'pressure_on','heat_on'};
% nCond       = 3;
% condName    = {'cue','redcross','whitecross'};
nPmod       = 0;%2;                %number parametric modulators
imageDim    = '4D';             %dimension of images to be processed (3D or 4D)

%first-level specific values
orthog          = 1;                    %as input for conditions (0 or 1)
skern           = 6;                    %for smooting
skernel         = repmat(skern,1,3);    %for smoothing
model           = 1;%5;%4;%3;%2;        %running number of model
modelName       = 'ramps';%'rampsUB';%'rampsCB';%'cues';%
basisF          = 'FIR'; %'HRF'         %one of HRF | FIR
units           = 'scans';              %onsets are specified in 'scans' or 'secs'
cvi             = 'none';               % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1)'|'AR(1) + w'
hpf             = 128;% 360;           % low frequency confound: high-pass cutoff (secs), Inf = no filtering, default is 128
nBase           = 12;%20;%              % depending on duration of longest stimulus (for fir model) 

% adjust the number of processors to work with to number of subjects------
if size(allSubs) < nProc
    nProc = size(allSubs,2);
end

% assign subjects to processors-------------------------------------------
%NOTE: splitSubs contain 1 x n_proc cells, each containing 1 or more subjects
splitSubs  = splitvect(allSubs, nProc);

%% define what batches to run
fprintf('\n##############################\nplease check your order now:\n##############################\n');
smoothie                = 1 %in case of smoothing realigned niftis: this has to be done first (separately)
disp('##############################');

% this has to be done separately
prepNormAndOrSmooth     = 0 %need this to do smooth and nl_smooth

%from UB-----------------------------------
norm2mni                = 0 %for now using those below from CB

%from CB-----------------------------------
% warp                    = 0 %using nl_warp
% smooth                  = 0 %using nl_smooth
nl_warp                 = 0
% smooth only works when done separately
nl_smooth               = 0


fprintf('\ncontinue with F5 or cancel with dbquit\n######################################\n');


%% loop across subjects ##################################################
% for np = 1:size(allSubs,2)
for np = 1:size(splitSubs,2)
    matlabbatch = [];
    
    % loop across subjects whitin each splitSubs-cell---------------------
    for g = 1:size(splitSubs{np},2)
        
        % define subjectspecific values-----------------------------------
        subName = sprintf('sub-%02.2d',splitSubs{np}(g));
        epis = struct('name',subName);
        hdrDir  = [derivDir 'spm_preprocessing' filesep subName filesep sprintf('ses-%02.2d',session) filesep];
        %-----------------------------------------------------------------
        fprintf(['\n\n SMOOTHING volunteer ' subName '\n']);
        
        
        % read in values from dicom header
        load(fullfile(hdrDir,[subName '_dcmhdr.mat']),'dcmhdr');
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
        firstLevelName  = sprintf('%s_%s%02d_%s_nBase%d_orth%d_tA%d',...
            studyName,modelName,model,basisF,nBase,orthog,adjValue);
        % firstLevelName = 'test';
        fprintf('\nbasis function is %s\norto is %d\ntempAdj is %d\n',basisF,orthog,tempAdj);
        
        %create empty variables-------------------------------------------
        epiDir     = cell(size(allRuns,2),1);
        epiFiles   = cell(size(epiDir,1),0);
        
        
        %looping across runs to get the files per run---------------------
        for nr = 1:size(allRuns,2)
            runName        = sprintf('run-%02.2d',nr);
            strnameRun     = sprintf('run_%02.2d',nr);
            epiDir         = fullfile(derivDir, 'spm_preprocessing',filesep, subName, sprintf('ses-%02.2d',session), subDir{2});
            %behDir         = fullfile(baseDir, derivDir, subName, subDir{5});
            anatDir        = fullfile(derivDir, 'spm_preprocessing',filesep, subName, sprintf('ses-%02.2d',session), subDir{3});
            
            %% Specify first level----------------------------------------
            cd(epiDir);%go to location of epi-files
            use_sNIFs = 0;
            if use_sNIFs == 1
                epiName = ['s' num2str(skern) '_ra' subName sprintf('_ses-%02.2d',session) '_task-' studyName '_' runName '_bold.nii'];
            else
                epiName = ['ra' subName sprintf('_ses-%02.2d',session) '_task-' studyName '_' runName '_bold.nii'];
            end
            
            numberScans = size(spm_select('expand',fullfile(epiDir,epiName)),1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %get image files
            %=========================================================================
            if strcmp(imageDim, '4D')
                epiFiles = cell(numberScans,1);
                
                for j = 1:numberScans
                    epiFiles{j}   = [epiDir filesep epiName ',' num2str(j)];
                end
                clear j
                
            elseif strcmp(imageDim, '3D') %@UB: not yet tested this version
                %                 epiFilter        = epiName;
                %                 epiDir           = [epiDir Files3D];
                %                 epiFiles         = char(select_scans(epiFilter, epiDir));
                
            end
            
            if isempty(epiFiles)
                fprintf('Attention: NO epi files!');
                return;
            end
            
            fprintf('\n%s has %d files in %s\n',subName,length(epiFiles),runName);
            epis.(strnameRun) = epiFiles;%each epi-name-string is one cell
            
           
            %get realignment files
            %=========================================================================
            rpFileName.(strnameRun) = fullfile(epiDir, ['rp_a' subName sprintf('_ses-%02.2d',session) '_task-' studyName '_' runName '_bold.txt']);
            
        end %for the loop across runs
        clear nr runName;
        
%         %setup epis for loop in matlabbatch
%         if size(allRuns,2) == 1
%             epis.run = {epis.run_01};
%             rpFileName.runs = {rpFileName.run_01};
%         else
%             epis.run = [{epis.run_01} {epis.run_02} {epis.run_03} {epis.run_04}];
%             rpFileName.runs = [{rpFileName.run_01} {rpFileName.run_02} {rpFileName.run_03} {rpFileName.run_04}];
%         end
        
%         %% set up 1level directory
%         %=========================================================================
%         %define name of firstlevel-dir
%         firstLevelDir = [derivDir 'spm_firstlevel' filesep subName filesep sprintf('ses-%02.2d',session) filesep];
%         
        %set index for matlabbatch
        gi   = 1;
        
        if smoothie == 1
            % smooth niftis BEFORE specifying firstlevel
            %=============================================================
            %get realigned niftis-----------------------------------------
            for nr = 1:size(allRuns,2) %@UB: the loop works when doing one sub at a time. I might remove this section from the firstlevel script
              
                
                runName     = sprintf('run-%02d',nr);
                funcDir     = fullfile(derivDir, 'spm_preprocessing', subName, sprintf('ses-%02.2d',session), subDir{2});
                ntfDir      = funcDir; %[funcDir filesep runName filesep];
                rNIFs   = char(spm_select('ExtFPList', ntfDir, ['^ra' subName sprintf('_ses-%02.2d',session) '_task-' studyName '_' runName '_bold.nii'])); % slicetimed and realigned niftis (rasub-01_task-pdev_run-01_bold.nii)              
           
                
                matlabbatch{gi}.spm.spatial.smooth.data   	= cellstr(rNIFs);
                matlabbatch{gi}.spm.spatial.smooth.fwhm   	= skernel;
                matlabbatch{gi}.spm.spatial.smooth.dtype  	= 0;
                matlabbatch{gi}.spm.spatial.smooth.im       = 0;
                matlabbatch{gi}.spm.spatial.smooth.prefix 	= ['s',num2str(skern) '_'];
                
                gi = gi + 1;
            end
        end %for smoothing

    end 
end 

if ~isempty(matlabbatch)
    spm_jobman('run',matlabbatch)
end


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

function [files] = select_scans(file_filt, direc)
x     = []; x     = spm_select('List', direc, file_filt);
y     = []; y     = [repmat([direc filesep], size(x, 1), 1) x repmat(',1', size(x, 1), 1)];
files = []; files = mat2cell(y, ones(size(y, 1), 1), size(y, 2));

function f_files = create_func_files(s_dir,f_templ,n_files)
for i=1:n_files
    f_files{i} = [s_dir filesep f_templ ',' num2str(i)];
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

        