function extract_participants_file
% BIDS format and generate json files and tsv files as well as events files

% based on scipt by Sebastian Puschmann (2018), University of Oldenburg
% adapted by Janne Nold (2022), UKE
%clear all

close all
subIDs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 
exclude_subs = [8];
session = 1;


if ~isempty(exclude_subs)
    subIDs = subIDs(~ismember(subIDs,exclude_subs));
end


% load in fitness and expectation file and choose correct subject
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "subject", "expectation_exercise","expectation_exercise_rating", "sporty","ftp","pwc"];
opts.SelectedVariableNames = ["subject", "expectation_exercise", "expectation_exercise_rating", "sporty","ftp","pwc"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
expectationsportydata = readtable("/projects/crunchie/nold/PEEP/behavioural/MAIN/expectation_sporty_data.csv", opts);
clear opts
exp_sporty= table2array(expectationsportydata);
exp_sporty = exp_sporty(2:end,:);


for subjects = 1:length(subIDs)

    %% Settings
    path        = get_study_specs;

     if ~exist([path.baseDir,'/rawdata/'])
        mkdir([path.baseDir,'/rawdata/'])
     end

    path2BIDS = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/';
    dataset_name = ['PEEP_MAIN'];
   

    %% Generate JSON and TSV on DATASET LEVEL
    % --------------------------------------------

    disp(['Generating BIDS JSON and TSV files at ' path2BIDS])

    % Generate dataset json
    % ----------------------------------------

    fid=fopen([path2BIDS filesep 'dataset_description.json'],'w');
    fprintf(fid,'{\n');
    fprintf(fid,['"Name": "' dataset_name '",\r\n']);
    fprintf(fid,['"BIDSVersion":  "1.1.1",\r\n']);
    fprintf(fid,['"Authors":  "Janne Nold",\r\n']);
    fprintf(fid,'}');
    fclose(fid);

    % Generate participants tsv file
    % ---------------------------------
    clear params
    cd(['/projects/crunchie/nold/PEEP/behavioural/MAIN' filesep sprintf('sub-%02.2d',subIDs(subjects)) filesep sprintf('ses-%02.2d',session) filesep]);
    params = load(sprintf('parameters_sub0%02.2d.mat',subIDs(subjects)));

    ind = find(exp_sporty(:,1) == subIDs(subjects));
    exp = exp_sporty(ind,2);
    exp_rating = exp_sporty(ind,3);
    fit = exp_sporty(ind,4);
    ftp_watt = exp_sporty(ind,5);
    pwc = exp_sporty(ind,6);

    if strcmp(params.P.subject.gender,'f')
        gender = 1;
    elseif strcmp(params.P.subject.gender,'m')
        gender = -1;
    end

    if ~exist([path2BIDS filesep 'participants.tsv'],'file')
        fid3=fopen([path2BIDS filesep 'participants.tsv'],'a+');
        fprintf(fid3,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s','participant_id','age','sex','fitness','ftp_watt','expectation','expecation_rating','relative_ftp');
        fprintf(fid3,'\r\n');
        fprintf(fid3,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f',sprintf('sub-%02.2d',params.P.protocol.subID),params.P.subject.age,gender, fit,round(ftp_watt,1), exp,round(exp_rating,2),round(pwc,2));
        fprintf(fid3,'\r\n');

    else
        fid3=fopen([path2BIDS filesep 'participants.tsv'],'a+');
        fprintf(fid3,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f',sprintf('sub-%02.2d',params.P.protocol.subID),params.P.subject.age,gender, fit,round(ftp_watt,1), exp,round(exp_rating,2),pwc);
        fprintf(fid3,'\r\n');
    end
    fclose(fid3);


 
        cd([path.baseDir,'/rawdata/']);
end
end
