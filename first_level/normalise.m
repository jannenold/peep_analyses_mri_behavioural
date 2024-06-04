function normalise(allSubs,model,modelname)

addpath(genpath('/home/nold/Desktop/PEEP/fMRI_new/'));

%% call hostname and define paths and number of processing units
hostName    = char(getHostName(java.net.InetAddress.getLocalHost));

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
        FFXDir           = [baseDir 'derivatives/spm_firstlevel/' basisF '_' modelname filesep];
        nProc            = 8; % 8 parallel processes on revelations
        subDir            = {'anat','func','fmap'};
    otherwise
        error('host not recognized');
end


%% define general studyspecific values / names

%subjects / runs / conditions
% allRuns     = 1:4;
% nSess       = 8;
% nPmod       = 0;%2;                %number parametric modulators
% imageDim    = '4D';             %dimension of images to be processed (3D or 4D)
% 
% %first-level specific values
% orthog          = 1;                    %as input for conditions (0 or 1)
% skern           = 3;                    %for smooting
% skernel         = repmat(skern,1,3);    %for smoothing
% %model           = 1;%5;%4;%3;%2;        %running number of model
% %modelName       = 'ramps';%'rampsUB';%'rampsCB';%'cues';%
% %basisF          = 'HRF'; %'HRF'         %one of HRF | FIR
% units           = 'scans';              %onsets are specified in 'scans' or 'secs'
% cvi             = 'none';               % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1)'|'AR(1) + w'
% hpf             = 128;% 360;           % low frequency confound: high-pass cutoff (secs), Inf = no filtering, default is 128
% nBase           = 12;%              % depending on duration of longest stimulus (for fir model) 

% adjust the number of processors to work with to number of subjects------
if size(allSubs) < nProc
    nProc = size(allSubs,2);
end

% assign subjects to processors-------------------------------------------
%NOTE: splitSubs contain 1 x n_proc cells, each containing 1 or more subjects
splitSubs  = splitvect(allSubs, nProc);

%% define what batches to run
fprintf('\n##############################\nplease check your order now:\n##############################\n');

check = 0  % 0 to run, 1 to not run chosen batch


disp('###############################');
% this has to be done separately
prepNormAndOrSmooth     = 1 %need this to do smooth and nl_smooth

%from UB-----------------------------------
norm2mni                = 1 %for now using those below from CB

%from CB-----------------------------------
nl_warp                 = 1



%% loop across subjects ##################################################
% for np = 1:size(allSubs,2)
for np = 1:size(splitSubs,2)
    matlabbatch = [];
    
    % loop across subjects whitin each splitSubs-cell---------------------
    for g = 1:size(splitSubs{np},2)
        
         % define subjectspecific values-----------------------------------
        subName = sprintf('sub-%02.2d',splitSubs{np}(g));
        epis = struct('name',subName);
        hdrDir  = [derivDir];
        firstLevelDir = [FFXDir subName filesep];
        epiDir         = fullfile(derivDir, subName, subDir{2});
        
        % read in values from dicom header
        load(fullfile(hdrDir,'dcmhdr.mat'),'dcmhdr');
        pixelSpacing    = dcmhdr.acqpar.PixelSpacing;
        sliceThickness  = dcmhdr.acqpar.SliceThickness;
        vox             = [pixelSpacing' sliceThickness];
        VboundingBox    = [NaN NaN NaN; NaN NaN NaN];

        %set index for matlabbatch
        gi   = 1;

        if prepNormAndOrSmooth == 1
            %setup for norm2mni / warping and smoothing
            %=================================================================
            %Get flowfields: images with the prefix "u_rc1s*"----------
            structDir   = fullfile(derivDir,subName , subDir{1});
           
            % Adjust for T1 taken ins ses-01 or ses-02
            file_ff     = ['u_rc1' subName '_T1w.nii'];
           
            if ~exist([structDir,filesep,file_ff])
                file_ff = ['u_rc1' subName sprintf('_ses-%02.2d',2) '_T1w.nii'];
            end 
            
            Flowfield   = char(spm_select('FPList', structDir, file_ff));


            % get betas, cons, mask-----------------------------------
            BETAs    = char(spm_select('FPList', firstLevelDir, '^beta.*nii$'));   % 1stLevel betas
            CONs     = char(spm_select('FPList', firstLevelDir, '^con.*nii$'));    % 1stLevel cons
            fCONs    = char(spm_select('FPList', firstLevelDir, '^ess_.*nii$'));  % 1stLevel F-cons
            %fCONs    = char(spm_select('FPList', firstLevelDir, '^spmF_.*nii$'));  % 1stLevel F-cons
            tCONs    = char(spm_select('FPList', firstLevelDir, '^spmT_.*nii$'));  % 1stLevel T-cons
            MASKs    = char(spm_select('FPList', firstLevelDir, 'mask*'));         % 1stLevel masks
            
        end %for prepNormAndOrSmooth
        
        if norm2mni == 1
            
            matlabbatch{gi}.spm.tools.dartel.mni_norm.template              = {''};
            matlabbatch{gi}.spm.tools.dartel.mni_norm.data.subj.flowfield   = cellstr(Flowfield);
            matlabbatch{gi}.spm.tools.dartel.mni_norm.data.subj.images      = cellstr(char([fCONs;CONs]));
            matlabbatch{gi}.spm.tools.dartel.mni_norm.vox        = vox;         % voxel dim taken from image header
            matlabbatch{gi}.spm.tools.dartel.mni_norm.bb         = VboundingBox;% if NaN = same size as scanned, isomorph voxel size is prefered
            matlabbatch{gi}.spm.tools.dartel.mni_norm.preserve   = 0;         % change for VBM?
            matlabbatch{gi}.spm.tools.dartel.mni_norm.fwhm       = [0 0 0];   % spatial smoothing default is [8 8 8];but here [0 0 0] because images unsmoothed
            
            gi = gi + 1;
        end %for norm2mni
        
  
        
        if nl_warp == 1
            %using nlin coreg + DARTEL
            %=============================================================
            matlabbatch{gi}.spm.util.defs.comp{1}.def = {[epiDir filesep 'y_epi_2_template.nii']}; 
            matlabbatch{gi}.spm.util.defs.out{1}.pull.fnames = cellstr(char([fCONs;CONs]));
            matlabbatch{gi}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
            matlabbatch{gi}.spm.util.defs.out{1}.pull.interp = 4;
            matlabbatch{gi}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{gi}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            matlabbatch{gi}.spm.util.defs.out{1}.pull.prefix = 'w_nlco_dartel';
            
            gi = gi + 1;
        end %for nl_warp
        
       


    end 
end 

if ~isempty(matlabbatch)
    spm_jobman('initcfg');
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
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