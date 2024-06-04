function overlay_movement_stimulus

% function check_realign(all_sub_ids);
% displays 1st images from all sessions in coreg to see whether everything is OK
all_sub_ids          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 

% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
TR = 1.8;
modality        = {'heat','pressure'};

for sub = 1:n_subs

    sub_id      = all_sub_ids(sub);
    all_runs_mov = [];
    cnt = 1;

    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            
            % Load Rp files of all runs and session
            all_runs_mov{cnt,1} = load(spm_select('FPList', [BIDS.dir,filesep,sprintf('sub-%02d',sub_id),filesep,'func'], sprintf('^rp_asub-%02d_ses-%02d_task-peep_run-%02d_bold.txt$',sub_id,ses,run)));

            
            % Load corresponding stimulus onset files
            dataDir = '/projects/crunchie/nold/PEEP/behavioural/MAIN';
            f_stub = [dataDir,filesep,sprintf('sub-%02.2d/ses-%02.2d/',sub_id,ses),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_',sub_id,run,ses+1)];
            [bl,tr, tim, event]= textread([f_stub 'all_events.tsv'], '%d %d %f %s', 'headerlines', 1, 'delimiter','\t');
            [sub, day, bl_2, tr_2, stim, int_VAS, int_phys, resp, resp_o, rating,	rt]= textread([f_stub 'behav.tsv'], '%d %d %d %d %s %d %f %d %f %d %f', 'headerlines', 1, 'delimiter','\t');

            all_event_start = find(contains(event,'start_plateau'));  

            for mod = 1:numel(modality)
                    mask = contains(stim,modality{mod});
                    RES{cnt,mod}    = (tim(all_event_start(mask))/TR);
            end

             cnt = cnt + 1;
        

        end
    end
    
    % Overlay Stimulus and onset files for each run

    for i = 1:(vars.nRuns*vars.nSess)
        subplot(ses*run,1,i);
        plot(all_runs_mov{i});
        legend('x','y','z')
        xline(RES{i,1},Color='red',LineWidth=2);
        xline(RES{i,2},Color='blue',LineWidth=2);
        title(sprintf('Sub-%02d',sub_id));
        x0=10;
        y0=10;
        width=900;
        height=1200;
        set(gcf,'position',[x0,y0,width,height]);
        legend('Pressure','Heat');
    end

    fprintf('Checking Subject %d \nPress enter in command window to continue\n',sub_id);
    input('');
end


           
