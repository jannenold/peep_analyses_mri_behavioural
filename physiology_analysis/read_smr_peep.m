%% import Spike2 data in .smr format and save in .mat format
% 220208 received from ma.wolf@uke.de
% 220208 adjusted by u.bromberg@uke.de

% For readability in Ledalab, the .mat file must be in a specific format.
% It must contain a single variable named data that is a struct. The struct
% must have fields named conductance, time, timeoff, and event. event must
% be a struct array with fields time, nid, name, userdata. time and nid
% are numeric; the name fields are character vectors, the userdata fields
% can be empty []

function read_smr_peep(sub,session)
%inarg
%sub:   digit between 1 and 99;
addpath(genpath('/home/nold/Documents/sigTOOL'));
path.base = ['/projects/crunchie/nold/PEEP/physiology/MAIN/'];
path.preproc = [path.base,'analysis/'];
rawext = sprintf('sub0%02d',sub);
path.raw = [path.base filesep 'original' filesep];
% addpath(genpath(path.toolbox));
%taskname = 'task-peep';
%suffix = 'recording-spike_stim';
expectedTR = 1.8;
mkdir([path.base sprintf('sub-%02.2d',sub) filesep sprintf('ses-%02.2d',session) filesep]);
save_path = [path.base sprintf('sub-%02.2d',sub) filesep sprintf('ses-%02.2d',session) filesep];
%% Automatic or Manual Selection of Files
fprintf('\n\n\n++++ Starting the SMR processing ++++\n\n');
if nargin < 1
    [files] = uigetfile([path.raw '*.smr'], 'Select Expression file(s)', 'MultiSelect', 'on');
    disp('Selected the following files:')
else
    subname = ['sub0' num2str(sub, '%02d');];
    day = sprintf('ses-%02d',session);
    filelist = dir([path.raw day filesep subname '*.smr']);
    files = {filelist.name};
    disp(['Found the following files for ' subname ':'])
end
disp(files)
disp(' ');
files = cellstr(files);

%% File Loop
for i = 1:length(files)
    thisFile = files{i};
    disp(['Processing ' thisFile ' (file ' num2str(i) ' of ' num2str(length(files)) ')...']);    
    sess = cellfun(@str2double, regexp(thisFile, '(?<=sess-)\d*', 'match'));
    runname = sprintf('run-%02d',sess);
   
    %% check if output file already exists
    % e.g. sub-01_task-pdev_run-01_recording-scr_stim.mat
    subbie = sprintf('sub-%02d',sub);
    if exist([save_path subbie '_' sprintf('ses-%02.2d',session) '.mat'], 'file')
        choice = input(['File ' thisFile ' has already been processed, press ''y'' to do the processing anyway or any other key to skip this file\n'], 's');
        if ~strcmp(choice, 'y')
            continue
        end
    end
    
    %% import kcl file
    kclfile = ImportSMR([path.raw day filesep thisFile], [path.preproc subname '/']);
    spd = load(kclfile, '-mat');
    
    %% read fieldnames
    fields = fieldnames(spd);
    nChannels = sum(startsWith(fields, 'head'));
    headFields = fields(startsWith(fields, 'head'));
    chanFields = fields(startsWith(fields, 'chan'));
    
    for j = 1:nChannels
        chanName = spd.(headFields{j}).title;
        
        if strcmp(spd.(headFields{j}).channeltype, 'Continuous Waveform')
            data.(chanName) = double(spd.(chanFields{j}).adc)*spd.(headFields{j}).adc.Scale;
            
            % identify marker precision
            mrkprec = spd.(chanFields{j}).tim(2)/size(spd.(chanFields{j}).adc,1);
            
            samplingInterval = double(spd.(headFields{j}).adc.SampleInterval(1))*double(spd.(headFields{j}).adc.SampleInterval(2));
            samplingRate = 1/(double(spd.(headFields{j}).adc.SampleInterval(1))*double(spd.(headFields{j}).adc.SampleInterval(2)));
            
        else  
            data.(chanName) = double(spd.(chanFields{j}).tim/mrkprec); 
        end
    end
    
    data.time = (1:length(data.SCR))'*samplingInterval;
    %data.pulse = data.PULS;
    %data.resp = data.RESP;
    data = rmfield(data, 'SCR');
    data.timeoff = 0;
    
    %% save timing of first scanner impulse
    data.firstMRPulse = data.scanner(1);
    
 
    %% Plot Repiration/Pulse and markers
    figure
    subplot(2,1,1)
    plot(data.Resp);
    xline(data.PainOn);
    xline(data.BlockOn,'LineStyle','-');
    %xline(data.RestOn);

    %xline(data.scanner);

    subplot(2,1,2)
    plot(data.PULS);
    xline(data.PainOn);
    xline(data.BlockOn,'LineStyle','-');
    %xline(data.RestOn);

    %xline(data.scanner);

%     for b = 1:4
%         data_new(b).PULS = data.PULS(data.BlockOn(b):data.BlockOn(b) + 600);
%         data_new(b).Resp = data.Resp(data.BlockOn(b):data.BlockOn(b) + 600);
%     end

    %% readLogFile for events
%     out = struct(readLogFile(sub,i));
% 
%    
%     %% setup structure 'data.event'
%     cnt = 1;   
%     for m = 1:size(out.log,2)
%        tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),1) = out.log{m}.onset';
%        switch out.log{m}.name
%            case 'flat-flat'
%                tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),2) = repmat(ones,size(out.log{m}.onset,2),1);
%            case 'steep-steep'
%                tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),2) = repmat(2,size(out.log{m}.onset,2),1);
%            case 'flat-steep'
%                tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),2) = repmat(3,size(out.log{m}.onset,2),1);
%            case 'steep-flat'
%                tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),2) = repmat(4,size(out.log{m}.onset,2),1);
%            case 'flat-flat_r'
%                tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),2) = repmat(5,size(out.log{m}.onset,2),1);
%            case 'steep-steep_r'
%                tmp(cnt:cnt+(size(out.log{m}.onset,2)-1),2) = repmat(6,size(out.log{m}.onset,2),1);
%        end
%        cnt = cnt + size(out.log{m}.onset,2);
%     end
%     clear out;
%     
%     tmpS = sortrows(tmp,1);
%     tmpD.time = tmpS(:,1)';
%     tmpD.nid = tmpS(:,2)';
%     
%     for n = 1:size(tmpD.time,2)
%        switch tmpD.nid(n)
%            case 1
%                tmpD.name{n} = 'flat-flat';
%            case 2
%                tmpD.name{n} = 'steep-steep';
%            case 3
%                tmpD.name{n} = 'flat-steep';
%            case 4
%                tmpD.name{n} = 'steep-flat';
%            case 5
%                tmpD.name{n} = 'flat-flat_r';
%            case 6
%                tmpD.name{n} = 'steep-steep_r';
%        end
%     end  
%     tmpD.userdata = [];
%     
%     for p = 1:size(tmpD.time,2)
%     data.event(p) = struct(...
%             'time',tmpD.time(p),...
%             'nid',tmpD.nid(p),...
%             'name',tmpD.name(p),...
%             'userdata',tmpD.userdata); 
%     end
%     
% %     %##############################################################################
% %         %% generate usable triggers
% %     
% %         % Cues:
% %         % 32: Cue Green
% %         % 64: Cue Red
% %         % 128: Cue Yellow
% %     
% %         % 96 - PreStim (all conditions)
% %         % 160 - PostStim (all conditions)
% %     
% %         % 192 - Start
% %         % 224 - End
% %     
% %         % 32 = TriggerA active
% %         % 64 = TriggerB active
% %         % 128 = TriggerC active
% %     
% %         tmp.tableA = table(data.TriggerA, repmat(32, length(data.TriggerA),1),...
% %             'VariableNames', {'time', 'TriggerA'});
% %         tmp.tableB = table(data.TriggerB, repmat(64, length(data.TriggerB),1),...
% %             'VariableNames', {'time', 'TriggerB'});
% %         tmp.tableC = table(data.TriggerC, repmat(128, length(data.TriggerC),1),...
% %             'VariableNames', {'time', 'TriggerC'});
% %         tmp.tableScanner = table(data.scanner, repmat(5, length(data.scanner),1),...
% %             'VariableNames', {'time', 'scanner'});
% %     
% %         triggerTable = outerjoin(...
% %             outerjoin(...
% %             outerjoin(...
% %             tmp.tableA, tmp.tableB,'MergeKeys', true),...
% %             tmp.tableC, 'MergeKeys', true),...
% %             tmp.tableScanner, 'MergeKeys', true);
% %     
% %         triggerTable.triggerValue = sum(...
% %             [...
% %             triggerTable.TriggerA,...
% %             triggerTable.TriggerB,...
% %             triggerTable.TriggerC,...
% %             triggerTable.scanner], 2, 'omitnan');
% %     
% %         clear tmp;
% %     
% %         eventLabels = {'Cue GREEN', 'Cue RED', 'Cue YELLOW', 'Start', 'End', 'scanner'};
% %         eventValues = [32, 64, 128, 192, 224, 5];
% %     
% %         for i = 1:length(eventLabels)
% %             triggerTable.Event(triggerTable.triggerValue == eventValues(i)) = eventLabels(i);
% %         end
% %     
% %         for i = 1:length(triggerTable.Event)
% %             if ismember(triggerTable.triggerValue(i), [32, 64, 128])
% %                 condition = regexp(triggerTable.Event{i}, '(?<=\s)\w+', 'match');
% %                 PreStimLine = find(triggerTable.triggerValue(i:end) == 96, 1) + i -1;
% %                 PainLine = find(triggerTable.triggerValue(i:end) == 160, 1) + i -1;
% %     
% %                 if ~isempty(PreStimLine)
% %                     triggerTable.Event{PreStimLine} = ['PreStim ' condition{1}];
% %                 end
% %                 if ~isempty(PainLine)
% %                     triggerTable.Event{PainLine} = ['Pain ' condition{1}];
% %                 end
% %             end
% %         end
% %     
% %         triggerTable.realTime = triggerTable.time/samplingRate;
% %         data.event = table2struct(...
% %             table(triggerTable.realTime, triggerTable.triggerValue, triggerTable.Event, double.empty(height(triggerTable),0), 'VariableNames', {'time', 'nid', 'name', 'userdata'}));
% %     
% %         %% check MR trigger timing
% %         mrtriggertimes = struct2table(data.event).time(strcmp(struct2table(data.event).name, 'scanner'));
% %         recordedTRs = diff(mrtriggertimes);
% %         totalDeviationFromTrueTR = sum(recordedTRs - expectedTR);
% %         if totalDeviationFromTrueTR > 0.1
% %             warning(['Recorded MR triggers for ' subname ' deviate strongly from the expected TR (' num2str(totalDeviationFromTrueTR) ' s). Please check!'])
% %         end
% %     

 

    %% save data
    suffix_new = [subbie '_' sprintf('ses-%02.2d',session) '.mat'];
    save([save_path '/' suffix_new], 'data');
    savefig([save_path '/' 'resp_pulse.fig']);
end

%#Subfunctions#############################################################
function [out] = readLogFile(sub,run)
%u.bromberg@uke.de
%sub: digit between 1:99
%run: digit between 1:3

[~, hostname] = system('hostname');
hostname = deblank(hostname);

if strcmp(hostname,'revelations')
    baseDir = '/projects/crunchie/nold/PEEP/fMRI/Data';
else
    baseDir = '/vicepa/imagen2/bromberg/cases/pdev/';
end

allRuns = 1:4;
subName = sprintf('sub-%02d',sub);
fprintf('\nrunning %s\n', subName);

if strcmp(hostname,'revelations')
    derivDir = [baseDir filesep studyPart filesep 'derivatives/'];
else
    derivDir = [baseDir 'data\' subName '\'];
end

cond_names = {'pressure_on','heat'};
out.log = [];
% out(numel(allRuns),numel(cond_names)).name = [];

% for j = 1:length(allRuns)
runName = sprintf('run-%02d',allRuns(run));
fprintf('\nrunning %s\n', runName);

%1. read in log-file
beh_file = ['^sub-\d\d_' runName '_task-pdev_log_.*\.mat'];
onset_file = spm_select('FPList',derivDir, beh_file);

%2. extract onsets
out.log = extract_offset(onset_file,cond_names);
% end

%---------------------------------------------------------------------------
function RES = extract_offset(onset_file,conditions)
i1=1;i2=1;i3=1;i4=1;
r1 = load(onset_file);
nl = size(r1.p.log.events,1);%number of events
for j=1:nl
    entry = r1.p.log.events(j,5);
    if iscell(entry{1}(2))
        switch entry{1}{2}
            case conditions{1}
                RES{1}.onset(i1) = r1.p.log.events{j,3};i1=i1+1; %read onset and count on
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
    RES{i}.name  = conditions{i}; %add name of condition
end
% finally take 1st and last ss ff as rating trials
RES{5}.onset = RES{1}.onset([1 end]);RES{1}.onset([1 end]) = [];RES{5}.name = 'flat-flat_r';
RES{6}.onset = RES{2}.onset([1 end]);RES{2}.onset([1 end]) = [];RES{6}.name = 'steep-steep_r';

%UB to-do##################################################################
%do I want the path. sub-fields?
%set up path subfunction for project to be used in all scripts






