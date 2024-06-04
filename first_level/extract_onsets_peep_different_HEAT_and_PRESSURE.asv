function RES = extract_onsets_peep_different_HEAT_and_PRESSURE(subjects,session,studyPart,nr,model)

modality        = {'heat','pressure'};
intensity       = {30,50,70};
TR              = 1.8;

dataDir = '/projects/crunchie/nold/PEEP/behavioural';
f_stub = [dataDir,filesep,studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_',subjects,nr,session+1)];

if nargin<1
    disp('No TSV file specivied. Using sub-01, ses-01');
    f_stub = '/projects/crunchie/nold/PEEP/behavioural/MAIN/sub-01/ses-01/sub-01_block-01-peep_exp-day-02_';
end

[bl,tr, tim, event]= textread([f_stub 'all_events.tsv'], '%d %d %f %s', 'headerlines', 1, 'delimiter','\t');
[sub, day, bl_2, tr_2, stim, int_VAS, int_phys, resp, resp_o, rating,rt]= textread([f_stub 'behav.tsv'], '%d %d %d %d %s %d %f %d %f %d %f', 'headerlines', 1, 'delimiter','\t');


%% Define different Onsets for Heat and Presure Stimuli
% Heat Onset is Plateau
% Pressure onset is stimulus onset (ca. Plateau - 1 TR)

start = {'start_plateau_heat','start_pressure'};
stop = {'stop_heat','stop_pressure'};

all_event_start = find(contains(event,start));   %start_plateau
all_event_stop  = find(contains(event,stop));    % stop_plateau for pressure and heat

cond = 1;

for mod = 1:numel(modality)
    for int = 1:numel(intensity)
        mask = contains(stim,modality{mod}) & (int_VAS == intensity{int});
        RES(cond).name     = [modality{mod} '_' num2str(intensity{int})];
        RES(cond).onset    = (tim(all_event_start(mask))/TR); 

        
        %RES(cond).duration = tim(all_event_stop(mask))/TR - tim(all_event_start(mask))/TR;
        %RES(cond).ratings  = rating(mask);
        %RES(cond).int_phys = int_phys(mask);
       
        if strcmp(model,'HRF')
            RES(cond).duration = tim(all_event_stop(mask))/TR - tim(all_event_start(mask))/TR;
            RES(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            RES(cond).tmod = 0;
            RES(cond).orth = 1;
        elseif strcmp(model,'FIR')
            RES(cond).duration = 0;
        end 

        cond = cond + 1;

end


end



end

