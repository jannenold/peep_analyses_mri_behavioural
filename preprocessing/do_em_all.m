function do_em_all


addpath('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/spm_bids_pipeline');

% define subjects that should be included
all_sub_ids          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 



do_import   = 0;
do_preproc  = 0;
do_qc       = 1;
do_analysis = 0;
do_means    = 0;
do_particpant_tsv = 0;

if do_import == 1
    import_data_ses;
end

if do_particpant_tsv == 1
    % if  aprticipants tsv file and json file should be generated in rawdata
    extract_participants_file
end

if do_preproc == 1
    % perform slice timing correction
    slice_timing(all_sub_ids);
    
    % realign part 1 of 2
    realign_1_2(all_sub_ids);
    
    % perform non-linear coregistartion to get from EPI space to T1 space
    nonlin_coreg(all_sub_ids);
    
    % perform skull stripping
    skullstrip(all_sub_ids);
    
    % perform spatial normalization of T1 to Template space using DARTEL
    create_dartel(all_sub_ids);
    
    % create various warp fields to map from EPI --> T1 --> Template space
    create_trans(all_sub_ids);
    
    % Create a masks for the 1st level GLM and brainstem
    create_mask(all_sub_ids);
    
    % realignemnt second step
    realign_2_2(all_sub_ids);
    
    % Warp a few images to Template space (skullstrip and mean EPI)
    warp_images(all_sub_ids);
        
    % Create 6 WM and 6 CSF noise regressors (and 6 regressors picking up noise
    % from the posterior part of the lateral ventricle)
    create_noise(all_sub_ids);
end

if do_means == 1
    % Create a mean skullstrip and a mean of mean EPIs
    create_means(all_sub_ids);
end

if do_qc == 1    
    % use coreg to display images for QC
    check_images(all_sub_ids);
    % Display realignment parameters (same as SPM) and check 1st image from all
    % runs
    check_realign(all_sub_ids);
end

if do_analysis == 1
    % do analyses
    analyses(all_sub_ids);    
end
end


