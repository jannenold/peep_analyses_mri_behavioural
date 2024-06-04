function RES = extract_onsets_peep_ONSET_START_STIMULUS_HEAT_and_PRESSURE(subjects,session,nr,model)

modality        = {'heat','pressure'};
intensity       = {30,50,70};
pharm_cond      = {'nacl','nax'};
TR              = 1.8;

% Pharma condition Load in
pharm_conds = '/home/nold/Desktop/PEEP/fMRI_new/scripts/first_level/PEEP_entblindungsliste_dummy_coded_final.csv';                
[subs,treatment_day2, treatment_day3, treatment_order]= textread(pharm_conds, '%d %d %d %s', 'headerlines', 1, 'delimiter',',');
treatment_order = str2num(cat(1,treatment_order{:}));

dataDir = '/projects/crunchie/nold/PEEP/behavioural';
f_stub = [dataDir,filesep,'MAIN',filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_',subjects,nr,session+1)];

[bl,tr, tim, event]= textread([f_stub 'all_events.tsv'], '%d %d %f %s', 'headerlines', 1, 'delimiter','\t');
[sub, day, bl_2, tr_2, stim, int_VAS, int_phys, resp, resp_o, rating,rt]= textread([f_stub 'behav.tsv'], '%d %d %d %d %s %d %f %d %f %d %f', 'headerlines', 1, 'delimiter','\t');


%% Define different Onsets for Heat and Presure Stimuli


start = {'start_heat','start_pressure'};
stop = {'stop_plateau_heat','stop_plateau_pressure'};

all_event_start = find(contains(event,start));   %start stimulus
all_event_stop  = find(contains(event,stop));    % stop_plateau for pressure and heat

startVAS = find(contains(event,'start_VAS'));
stopVAS = find(contains(event,'stop_VAS'));

cond = 1;
       


for mod = 1:numel(modality)
    for int = 1:numel(intensity)

        mask = contains(stim,modality{mod}) & (int_VAS == intensity{int});

        if treatment_order(find(subs == subjects)) == 0 && session == 1
            RES(cond).name     = [modality{mod} '_' num2str(intensity{int}) '_nacl'];
        elseif treatment_order(find(subs == subjects)) == 0 && session == 2
            RES(cond).name     = [modality{mod} '_' num2str(intensity{int}) '_nax'];
        elseif treatment_order(find(subs == subjects)) == 1 && session == 1
            RES(cond).name     = [modality{mod} '_' num2str(intensity{int}) '_nax'];
        elseif treatment_order(find(subs == subjects)) == 1 && session == 2
            RES(cond).name     = [modality{mod} '_' num2str(intensity{int}) '_nacl'];
        end

        RES(cond).onset    = (tim(all_event_start(mask))/TR);


        if strcmp(model,'HRF')
            RES(cond).duration = (tim(all_event_stop(mask))/TR - tim(all_event_start(mask))/TR);
            RES(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            RES(cond).tmod = 0;
            RES(cond).orth = 1;
        elseif strcmp(model,'FIR')
            RES(cond).duration = 0;
        end

        cond = cond + 1;
    end
end

if strcmp(model,'HRF')
    % Log VAS Rating as condition
    RES(cond).name  = 'VAS';
    RES(cond).onset    = (tim(startVAS)/TR);
    RES(cond).duration = tim(stopVAS)/TR - tim(startVAS)/TR;
    RES(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
    RES(cond).tmod = 0;
    RES(cond).orth = 1;
end
end



