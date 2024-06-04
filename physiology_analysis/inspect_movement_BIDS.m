function inspect_movement_BIDS(subs,nRun,nSess,insMov,plotOn,movReg,threshSpike,threshSess)

fprintf('\n\n--------------------------------------------------------\n');
fprintf('                       Spikes                               \n');
fprintf('------------------------------------------------------------\n');

    basedir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_preprocessing';
    outdir = '/projects/crunchie/nold/PEEP/physiology/MAIN/noise_regressors';
    percSpike = 0.05;

    if insMov == 1
    for g=1:numel(subs)
        subdir = fullfile(basedir,sprintf('sub-%02.0f',subs(g)),'func');
        mParams = [];
        for s = 1:nSess
            for j = 1:nRun
                if nSess == 1
                    rpFile = fullfile(subdir,sprintf('rp_asub-%02.0f_task-peep_run-%02.0f_bold.txt',subs(g),j));
                    matFile = fullfile(subdir,sprintf('asub-%02.0f_task-peep_run-0%d_bold.mat',subs(g),j));
                else
                    rpFile = fullfile(subdir,sprintf('rp_asub-%02.0f_ses-%02.0f_task-peep_run-%02.0f_bold.txt',subs(g),s,j));
                    matFile = fullfile(subdir,sprintf('asub-%02.0f_ses-%02.0f_task-peep_run-%02.0f_bold.mat',subs(g),s,j));
                end
                run_data = load(rpFile);
                load(matFile);
                maxMov = max(abs(run_data(:,1:3)))';
                if any(maxMov > threshSess)
                    maxM = max(maxMov(maxMov > threshSess));
                    fprintf('Sub%02.0f Run%d: movement within session of %.2f mm detected, please check participant\n',subs(g),j,maxM);
                end
                matFirst(:,:,nRun*s-(nRun-j)) = mat(:,:,2);
                matLast(:,:,nRun*s-(nRun-j)) = mat(:,:,end);
                withinSessMov = matLast(:,:,nRun*s-(nRun-j))/matFirst(:,:,nRun*s-(nRun-j));
                wsm(:,nRun*s-(nRun-j)) = abs(withinSessMov(1:3,4));
                if any(wsm > threshSess)
                    fprintf('Sub%02.0f Run%d: movement between first and last image of %.2f mm detected, please check participant\n',subs(g),j,max(wsm));
                end
                mParams = [mParams; run_data];
            end
        end
        
         for j = 1:nRun*nSess-1
                betSessMov(:,:,j) = matLast(:,:,j)/matFirst(:,:,j+1);
                mov2First(:,:,j) = matLast(:,:,j)/matFirst(:,:,1);
                bsm(:,j) = abs(betSessMov(1:3,4,j));
         end
            
        if any(bsm(:) > threshSess)
            fprintf('Sub%02.0f: movement between sessions of %.2f mm detected, please check participant\n',subs(g),max(bsm(:)));
        end
        if plotOn == 1
            figure;
            subplot(3,1,1)
            plot(mParams(:,1:3));legend('x','y','z');title(sprintf('Sub%02.0f',subs(g)));
            title(sprintf('Sub%02.0f',subs(g)));
            subplot(3,1,2)
            plot([zeros(3,1) bsm]','o-');legend('x','y','z');xlim([0.5 nRun*nSess+0.5]);ylim([-3 3]);title('Between session movement');
            subplot(3,1,3)
            plot(wsm','o-');xlim([0.5 nRun*nSess+0.5]);ylim([-3 3]);title('Within session movement (last-first)')
            d = input('Continue with next subject?');
        end
    end
    end
    
    if movReg == 1
        % find spike movement and create noise regressors for first level
        for g=1:numel(subs)
            subdir = fullfile(basedir,sprintf('sub-%02.0f',subs(g)),'func');
            for s = 1:nSess
                for j = 1:nRun
                    if nSess == 1
                        rpFile = fullfile(subdir,sprintf('rp_asub-%02.0f_task-peep_run-%02.0f_bold.txt',subs(g),j));
                        outFilename = sprintf('sub-%02.0f_nui_reg_mov_run%02.0f.mat',subs(g),j);
                    else
                        rpFile = fullfile(subdir,sprintf('rp_asub-%02.0f_ses-%02.0f_task-peep_run-%02.0f_bold.txt',subs(g),s,j));
                        outFilename = sprintf('sub-%02.0f_nui_reg_mov_ses-%02.0f_run-%02.0f.mat',subs(g),s,j);
                    end
                    mParams = load(rpFile);
                    diffParams = abs(diff(mParams(:,1:3)));
                    [nSpikes,c] = find(diffParams>threshSpike);
                    nSpikes = unique(nSpikes);
                    if numel(nSpikes) ~= 0
                        if numel(nSpikes) > floor(size(mParams,1)*percSpike)
                            fprintf('Sub%02.0f Run%d: spike movement in more than %.1f%% of volumes detected, please check participant\n',subs(g),j,percSpike*100);
                        end
                        fprintf('Sub%02.0f, Run%d: %d spikes detected\n',subs(g),j,numel(nSpikes));
                        nui_reg = zeros(size(mParams,1),numel(nSpikes));
                        for m = 1:length(nSpikes)
                            nui_reg(nSpikes(m)+1,m) = 1;
                        end
                        fprintf('Saving noise regressors for respective images\n');
                        save(fullfile(outdir,outFilename),'nui_reg');
                    end
                end
            end
        end
    end
    
   