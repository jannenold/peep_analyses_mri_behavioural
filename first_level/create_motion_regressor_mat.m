function create_motion_regressor_mat(subs,session)
clear R
blocks = 1:4;


for b = 1:length(blocks)
 fid = fopen(['/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/derivatives/spm_preprocessing/' sprintf('sub-%02.2d',subs) filesep sprintf('ses-%02.2d',session) '/func/',sprintf('rp_asub-%02.2d_ses-%02.2d_task-peep_run-%02.2d_bold.txt',subIDs,session,blocks(b))]);
    R2{b,1} = textscan(fid,'%f %f %f %f %f %f');
    fclose(fid);  
end 

R = [];
for i = 1:length(blocks)
R(i,:) = cell2mat(R2{i});
end 

%save(['/projects/crunchie/nold/PEEP/fMRI/Data/' studyPart '/derivatives/spm_firstlevel/' sprintf('sub-%02.2d',subIDs) filesep sprintf('ses-%02.2d',session) filesep 'R.mat'],'R');