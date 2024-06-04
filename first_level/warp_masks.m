function warp_masks(subjects,model,modelname)
%function warp_masks
%warp firstlevel masks to common space. second level mask is based on this
% adapted from Lukas Neugebauer, 
addpath(genpath('/home/nold/Desktop/PEEP/fMRI/'));

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


n_subs = length(subjects);

matlabbatch = cell(n_subs, 1);
counter = 0;

for s = 1:length(subjects)

    sub = subjects(s);
    counter = counter + 1;
    fl_dir_sub  = fullfile(FFXDir, sprintf('sub-%02d',sub));
    anatDir    = fullfile(derivDir,  sprintf('sub-%02d',sub), 'anat');
    funcDir    = fullfile(derivDir,  sprintf('sub-%02d',sub), 'func');
    mask        = fullfile(fl_dir_sub, 'mask.nii');
    new_name    = ins_letter(mask, sprintf('sub-%02d_', sub));
    flowfield_epi = char(spm_select('FPList',funcDir,'y_epi_2_template.nii'));

    if ~exist(new_name, 'file')
        try
        movefile(mask, new_name);
        catch
            keyboard
        end
    end
    
    matlabbatch{counter}.spm.util.defs.comp{1}.def = {flowfield_epi};
    matlabbatch{counter}.spm.util.defs.out{1}.pull.fnames = cellstr(new_name);
    matlabbatch{counter}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{counter}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{counter}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{counter}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{counter}.spm.util.defs.out{1}.pull.prefix = 'w';


end

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);


% %loop again, move files
% for s = 1:n_subs
% 
%     sub = subjects(s);
%     fl_dir_sub = fullfile(FFXDir,  sprintf('sub-%02d',sub));
%     anatDir = fullfile(derivDir,  sprintf('sub-%02d',sub), 'anat');
%     funcDir = fullfile(derivDir,  sprintf('sub-%02d',sub), 'func');
%     mod_dir = fullfile(fl_dir_sub);
%     mask = fullfile(funcDir, sprintf('wsub-%02d_mask.nii', sub));
%     if ~exist(mask, 'file')
%         error('No mask for %d ', sub);
%     end
%     movefile(mask, mod_dir);
% 
% end

function out = ins_letter(pscan,letter_start,letter_end)
if nargin <3
    letter_end = [];
end
for a=1:size(pscan,1)
    [p , f, e] = fileparts(pscan(a,:));
    out(a,:) = [p filesep letter_start f letter_end e];
end
end
end
