function [onsets] = extract_onsetfile_peep(subjects,session,studyPart)

% Create Onsets in SCANS by dividing onsets in seconds/TR

blocks = 1:4;
TR = 1.8; % in seconds 
  

% Preallocate
    clear onsets;

    onsets(b,1).name = 'pressure';
    onsets(b,2).name = 'heat';
%   onsets(b,1).onset = [];
%    onsets(b,2).onset = [];
%    onsets(b,1).duration = 0;
%    onsets(b,2).duration = 0;


for b = 1:length(blocks)
    clear C
   
    fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
    C = textscan(fid,'%f %f %f %s','HeaderLines',1);
    fclose(fid);

  
for i = 1:size(C{1,1},1)

    if strcmp(C{1,4}{i,1},'start_pressure')

        onsets(b,1).onset = [onsets(b,1).onset;C{1,3}(i)/TR] ;

    elseif strcmp(C{1,4}{i,1},'start_heat')
      
        onsets(b,2).onset = [onsets(b,2).onset;C{1,3}(i)/TR] ;

    end 
end

end 


mkdir(['/projects/crunchie/nold/PEEP/fMRI/Data/' studyPart '/derivatives/spm_firstlevel/' sprintf('sub-%02.2d',subjects) filesep sprintf('ses-%02.2d',session)]);
save(['/projects/crunchie/nold/PEEP/fMRI/Data/' studyPart '/derivatives/spm_firstlevel/' sprintf('sub-%02.2d',subjects) filesep sprintf('ses-%02.2d',session) filesep 'onsets.mat'],'onsets')

