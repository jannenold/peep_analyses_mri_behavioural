function compute_2ndlevel_mask(model,subs,exclude_subs,modelname)
%
%Compute second level mask on basis of mask.nii from firstlevel masks

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
        RFXDir           = [baseDir 'derivatives/spm_secondlevel']; % untergeordnet
       

    otherwise
        error('host not recognized');
end

if ~exist(RFXDir,"dir")
mkdir(RFXDir);
end



if ~isempty(exclude_subs)
    subs = subs(~ismember(subs,exclude_subs));
end

n_subs    = length(subs);


p_valid = 1; % 80% of the subs should have a 1 there



mask_templ = '^w.*mask\.nii$';
masks = cell(length(subs), 1);

for s = 1:n_subs

    sub         = sprintf('sub-%02d',subs(s));
    subdir      = fullfile(FFXDir, sub);
    m           = spm_select('FPList', subdir, mask_templ);
    assert(size(m, 1) == 1);
    masks{s}    = m;

end

matlabbatch = {};
matlabbatch{1}.spm.util.imcalc.input          = masks;
matlabbatch{1}.spm.util.imcalc.output         = 'mask_secondlevel.nii';
matlabbatch{1}.spm.util.imcalc.outdir         = {RFXDir};
matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
% only take voxels where we have data from all subjects
matlabbatch{1}.spm.util.imcalc.expression     = 'all(X)';%'all(X)';%;sprintf('sum(X) > %d', floor(n_subs * p_valid));
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 1;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;


spm_jobman('run', matlabbatch);

end
