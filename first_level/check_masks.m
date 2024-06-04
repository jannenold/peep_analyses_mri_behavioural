function check_masks(all_sub_ids,model,modelname)

% function check_images(all_sub_ids);
% displays resulting images in coreg to see whether everything is OK


% resolve paths
template_path = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/spm_bids_pipeline/templates/';
n_subs        = length(all_sub_ids);

if model == 1
    basisF          = 'HRF'; %'HRF'         %one of HRF | FIR
elseif model == 2
    basisF          = 'FIR'; %'HRF'         %one of HRF | FIR
end

baseDir          = ['/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/']; % übergeordnet
FFXDir           = [baseDir 'derivatives/spm_firstlevel/' basisF '_' modelname filesep];



for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
   
    ffx_dir   = fullfile(FFXDir, sprintf('sub-%02d',sub_id));
    
    mask_file      = fullfile(ffx_dir, sprintf('wsub-%02d_mask.nii',sub_id));
    template_file   = fullfile(template_path,'cb_Template_T1.nii');
    
    matlabbatch{1}.spm.util.checkreg.data = cellstr([{mask_file};{template_file}]);
    % Images in the same row are in the same space !
    
    spm_jobman('run',matlabbatch);
    fprintf('Checking Subject %d \nPress enter in command window to continue\n',sub_id);
    input('');
end

end