function [onsets] = extract_onsets_peep(subjects,session,studyPart,nr,model,nConds)

% Create Onsets in SCANS by dividing onsets in seconds/TR
if strcmp(model,'FIR')

        blocks = nr;
        
        for b = 1:length(blocks)
            clear fid
            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_behav.tsv',subjects,blocks(b),session+1)]);
            D = textscan(fid,'%f %f %f %f %s %f %f %f %f %f %f','HeaderLines',1);

            % Extract trials for different pain intensity 
            int_70 = D{1,4}(D{1,6} == 70);
            int_50 = D{1,4}(D{1,6} == 50);
            int_30 = D{1,4}(D{1,6} == 30);
            
            fclose(fid);


        end

        TR = 1.8; % in seconds

        for b = 1:length(blocks)
            clear C

            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
            C = textscan(fid,'%f %f %f %s','HeaderLines',1);
            fclose(fid);

            % Preallocate
            clear onsets;

            onsets(b,1).name = 'pressure_70';
            onsets(b,2).name = 'pressure_50';
            onsets(b,3).name = 'pressure_30';
            onsets(b,4).name = 'heat_70';
            onsets(b,5).name = 'heat_50';
            onsets(b,6).name = 'heat_30';
            onsets(b,1).onset = [];
            onsets(b,2).onset = [];
            onsets(b,3).onset = [];
            onsets(b,4).onset = [];
            onsets(b,5).onset = [];
            onsets(b,6).onset = [];
            onsets(b,1).duration = 0;
            onsets(b,2).duration = 0;
            onsets(b,3).duration = 0;
            onsets(b,4).duration = 0;
            onsets(b,5).duration = 0;
            onsets(b,6).duration = 0;

  
            counter_70 = 1;
            counter_50 = 1;
            counter_30 = 1;

            for i = 1:size(C{1,1},1)

                % extract pressure  onsets 
                if strcmp(C{1,4}{i,1},'start_plateau_pressure') 
    
                    % extract 70 VAS onsets
                    if C{1,2}(i) == int_70(counter_70)
                       
                        onsets(b,1).onset = [onsets(b,1).onset;C{1,3}(i)/TR] ;
                        counter_70 = counter_70 +1;
                        if counter_70 == 7
                           counter_70 = 6;
                        end 

                    % extract pressure 50 VAS onset
                    elseif C{1,2}(i) == int_50(counter_50)
                     
                        onsets(b,2).onset = [onsets(b,2).onset;C{1,3}(i)/TR] ;
                        counter_50 = counter_50 +1;
                        if counter_50 == 7
                            counter_50 = 6;
                        end
                        
                        % extract pressure 30 VAS onset 
                    elseif C{1,2}(i) == int_30(counter_30)
                     
                        onsets(b,3).onset = [onsets(b,3).onset;C{1,3}(i)/TR] ;
                        counter_30 = counter_30 +1;
                        if counter_30 == 7
                            counter_30 = 6;

                        end
                    end


                elseif strcmp(C{1,4}{i,1},'start_plateau_heat')

                   % extract 70 VAS onsets
                    if C{1,2}(i) == int_70(counter_70)
                       
                        onsets(b,4).onset = [onsets(b,4).onset;(C{1,3}(i)/TR)] ;
                        counter_70 = counter_70 +1;
                        if counter_70 == 7
                            counter_70 = 6;
                        end 

                    % extract  50 VAS onset
                    elseif C{1,2}(i) == int_50(counter_50)
                     
                        onsets(b,5).onset = [onsets(b,5).onset;(C{1,3}(i)/TR)] ;
                        counter_50 = counter_50 +1;
                        if counter_50 == 7
                            counter_50 = 6;
                        end
                        % extract  30 VAS onset 
                    elseif C{1,2}(i) == int_30(counter_30)
                     
                        onsets(b,6).onset = [onsets(b,6).onset;(C{1,3}(i)/TR)] ;
                        counter_30 = counter_30 +1;
                        if counter_30 == 7
                            counter_30 = 6;

                        end
                    end
 
                end
            end

        end

        %save('onsets.mat',"onsets");


    elseif nConds <= 2
        blocks = nr;
        TR = 1.8; % in seconds

        for b = 1:length(blocks)
            clear C

            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
            C = textscan(fid,'%f %f %f %s','HeaderLines',1);
            fclose(fid);

            % Preallocate
            clear onsets;

            onsets(b,1).name = 'pressure';
            onsets(b,2).name = 'heat';
            onsets(b,1).onset = [];
            onsets(b,2).onset = [];
            onsets(b,1).duration = 0;
            onsets(b,2).duration = 0;



            for i = 1:size(C{1,1},1)

                if strcmp(C{1,4}{i,1},'start_pressure')

                    onsets(b,1).onset = [onsets(b,1).onset;C{1,3}(i)/TR] ;
                   
                elseif strcmp(C{1,4}{i,1},'start_heat')

                    onsets(b,2).onset = [onsets(b,2).onset;C{1,3}(i)/TR] ;
                   
                end
            end

        end

    


elseif strcmp(model,'HRF')

    if nConds <= 2
        blocks = nr;
        TR = 1.8; % in seconds

        for b = 1:length(blocks)
            clear C

            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
            C = textscan(fid,'%f %f %f %s','HeaderLines',1);
            fclose(fid);

            % Preallocate
            clear onsets;

            onsets(b,1).name = 'pressure';
            onsets(b,2).name = 'heat';
            onsets(b,1).onset = [];
            onsets(b,2).onset = [];
            onsets(b,1).duration = [];
            onsets(b,2).duration = [];
            onsets(b,1).tmod = 0;
            onsets(b,2).tmod = 0;
            onsets(b,1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,1).orth = 1;
            onsets(b,2).orth = 1;


            for i = 1:size(C{1,1},1)

                if strcmp(C{1,4}{i,1},'start_pressure')

                    onsets(b,1).onset = [onsets(b,1).onset;C{1,3}(i)/TR] ;
                    onsets(b,1).duration = [onsets(b,1).duration;((C{1,3}(i+1)/TR)-(C{1,3}(i)/TR))]; %subtracting stop plateau from start plateau

                elseif strcmp(C{1,4}{i,1},'start_heat')

                    onsets(b,2).onset = [onsets(b,2).onset;C{1,3}(i)/TR] ;
                    onsets(b,2).duration = [onsets(b,2).duration;((C{1,3}(i+1)/TR)-(C{1,3}(i)/TR))]; %subtracting stop plateau from start plateau

                end
            end

        end

    elseif nConds == 6
        
        % Extract onsets for each pain intensity and modality
        blocks = nr;
        
        for b = 1:length(blocks)
            clear fid
            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_behav.tsv',subjects,blocks(b),session+1)]);
            D = textscan(fid,'%f %f %f %f %s %f %f %f %f %f %f','HeaderLines',1);

%             % Pressure 70
%             idx_int_70 = D{1,4}(D{1,6} == 70);
%             idx_pressure = D{1,4}(strcmp(D{1,5},'pressure'));
%             pressure_70 = intersect(idx_int_70,idx_pressure);
%             
%            % pressure 50
%             idx_int_50 = D{1,4}(D{1,6} == 50);
%             idx_pressure = D{1,4}(strcmp(D{1,5},'pressure'));
%             pressure_50 = intersect(idx_int_50,idx_pressure);
%             
%             % Presusre 30
%             idx_int_30 = D{1,4}(D{1,6} == 30);
%             idx_pressure = D{1,4}(strcmp(D{1,5},'pressure'));
%             pressure_30 = intersect(idx_int_30,idx_pressure);
%             
%             % Heat 70
%             idx_int_70 = D{1,4}(D{1,6} == 70);
%             idx_heat = D{1,4}(strcmp(D{1,5},'heat'));
%             heat_70 = intersect(idx_int_70,idx_heat);
%             
%            % pressure 50
%             idx_int_50 = D{1,4}(D{1,6} == 50);
%             idx_heat = D{1,4}(strcmp(D{1,5},'heat'));
%             heat_50 = intersect(idx_int_50,idx_heat);
%             
%             % H 30
%             idx_int_30 = D{1,4}(D{1,6} == 30);
%             idx_heat = D{1,4}(strcmp(D{1,5},'heat'));
%             heat_30 = intersect(idx_int_30,idx_heat);
            
            % Extract trials for different pain intensity 
            int_70 = D{1,4}(D{1,6} == 70);
            int_50 = D{1,4}(D{1,6} == 50);
            int_30 = D{1,4}(D{1,6} == 30);
            
            fclose(fid);


        end


        
        blocks = nr;
        TR = 1.8; % in seconds

        for b = 1:length(blocks)
            clear C

            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
            C = textscan(fid,'%f %f %f %s','HeaderLines',1);
            fclose(fid);

            % Preallocate
            clear onsets;

            onsets(b,1).name = 'pressure_70';
            onsets(b,2).name = 'pressure_50';
            onsets(b,3).name = 'pressure_30';
            onsets(b,4).name = 'heat_70';
            onsets(b,5).name = 'heat_50';
            onsets(b,6).name = 'heat_30';
            onsets(b,1).onset = [];
            onsets(b,2).onset = [];
            onsets(b,3).onset = [];
            onsets(b,4).onset = [];
            onsets(b,5).onset = [];
            onsets(b,6).onset = [];
            onsets(b,1).duration = [];
            onsets(b,2).duration = [];
            onsets(b,3).duration = [];
            onsets(b,4).duration = [];
            onsets(b,5).duration = [];
            onsets(b,6).duration = [];
            onsets(b,1).tmod = 0;
            onsets(b,2).tmod = 0;
            onsets(b,3).tmod = 0;
            onsets(b,4).tmod = 0;
            onsets(b,5).tmod = 0;
            onsets(b,6).tmod = 0;
            onsets(b,1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,3).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,4).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,5).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,6).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,1).orth = 1;
            onsets(b,2).orth = 1;
            onsets(b,3).orth = 1;
            onsets(b,4).orth = 1;
            onsets(b,5).orth = 1;
            onsets(b,6).orth = 1;

  
            counter_70 = 1;
            counter_50 = 1;
            counter_30 = 1;

            for i = 1:size(C{1,1},1)

                % extract pressure  onsets 
                if strcmp(C{1,4}{i,1},'start_pressure') 
    
                    % extract 70 VAS onsets
                    if C{1,2}(i) == int_70(counter_70)
                       
                        onsets(b,1).onset = [onsets(b,1).onset;C{1,3}(i)/TR] ;
                        onsets(b,1).duration = [onsets(b,1).duration;((C{1,3}(i+3)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_70 = counter_70 +1;
                        if counter_70 == 7
                           counter_70 = 6;
                        end 

                    % extract pressure 50 VAS onset
                    elseif C{1,2}(i) == int_50(counter_50)
                     
                        onsets(b,2).onset = [onsets(b,2).onset;C{1,3}(i)/TR] ;
                        onsets(b,2).duration = [onsets(b,2).duration;((C{1,3}(i+3)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_50 = counter_50 +1;
                        if counter_50 == 7
                            counter_50 = 6;
                        end
                        
                        % extract pressure 30 VAS onset 
                    elseif C{1,2}(i) == int_30(counter_30)
                     
                        onsets(b,3).onset = [onsets(b,3).onset;C{1,3}(i)/TR] ;
                        onsets(b,3).duration = [onsets(b,3).duration;((C{1,3}(i+3)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_30 = counter_30 +1;
                        if counter_30 == 7
                            counter_30 = 6;

                        end
                    end


                elseif strcmp(C{1,4}{i,1},'start_heat')

                   % extract 70 VAS onsets
                    if C{1,2}(i) == int_70(counter_70)
                       
                        onsets(b,4).onset = [onsets(b,4).onset;C{1,3}(i)/TR] ;
                        onsets(b,4).duration = [onsets(b,4).duration;((C{1,3}(i+3)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_70 = counter_70 +1;
                        if counter_70 == 7
                            counter_70 = 6;
                        end 

                    % extract pressure 50 VAS onset
                    elseif C{1,2}(i) == int_50(counter_50)
                     
                        onsets(b,5).onset = [onsets(b,5).onset;C{1,3}(i)/TR] ;
                        onsets(b,5).duration = [onsets(b,5).duration;((C{1,3}(i+3)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_50 = counter_50 +1;
                        if counter_50 == 7
                            counter_50 = 6;
                        end
                        % extract pressure 30 VAS onset 
                    elseif C{1,2}(i) == int_30(counter_30)
                     
                        onsets(b,6).onset = [onsets(b,6).onset;C{1,3}(i)/TR] ;
                        onsets(b,6).duration = [onsets(b,6).duration;((C{1,3}(i+3)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_30 = counter_30 +1;
                        if counter_30 == 7
                            counter_30 = 6;

                        end
                    end
 
                end
            end

        end

    elseif nConds == 3

        % Extract onsets for each pain intensity
        blocks = nr;
        
        for b = 1:length(blocks)
            clear fid
            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_behav.tsv',subjects,blocks(b),session+1)]);
            D = textscan(fid,'%f %f %f %f %s %f %f %f %f %f %f','HeaderLines',1);

            % Extract trials for different pain intensity 
            int_70 = D{1,4}(D{1,6} == 70);
            int_50 = D{1,4}(D{1,6} == 50);
            int_30 = D{1,4}(D{1,6} == 30);


            fclose(fid);


        end

        
        blocks = nr;
        TR = 1.8; % in seconds

        for b = 1:length(blocks)
            clear C

            fid = fopen(['/projects/crunchie/nold/PEEP/behavioural/',studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
            C = textscan(fid,'%f %f %f %s','HeaderLines',1);
            fclose(fid);

            % Preallocate
            clear onsets;

            onsets(b,1).name = 'int_70';
            onsets(b,2).name = 'int_50';
            onsets(b,3).name = 'int_30';
            onsets(b,1).onset = [];
            onsets(b,2).onset = [];
            onsets(b,3).onset = [];
            onsets(b,1).duration = [];
            onsets(b,2).duration = [];
            onsets(b,3).duration = [];
            onsets(b,1).tmod = 0;
            onsets(b,2).tmod = 0;
            onsets(b,3).tmod = 0;
            onsets(b,1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,3).pmod = struct('name', {}, 'param', {}, 'poly', {});
            onsets(b,1).orth = 1;
            onsets(b,2).orth = 1;
            onsets(b,3).orth = 1;


            counter_70 = 1;
            counter_50 = 1;
            counter_30 = 1;


            for i = 1:size(C{1,1},1)

                if strcmp(C{1,4}{i,1},'start_pressure') || strcmp(C{1,4}{i,1},'start_heat') 
                    
                    % extract  70 VAS onsets
                    if C{1,2}(i) == int_70(counter_70)
                       
                        onsets(b,1).onset = [onsets(b,1).onset;C{1,3}(i)/TR] ;
                        onsets(b,1).duration = [onsets(b,1).duration;((C{1,3}(i+1)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_70 = counter_70 +1;
                        if counter_70 == 7
                            counter_70 = 6;
                        end 

                    % extract pressure 50 VAS onset
                    elseif C{1,2}(i) == int_50(counter_50)
                     
                        onsets(b,2).onset = [onsets(b,2).onset;C{1,3}(i)/TR] ;
                        onsets(b,2).duration = [onsets(b,2).duration;((C{1,3}(i+1)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_50 = counter_50 +1;
                        if counter_50 == 7
                            counter_50 = 6;
                        end
                        % extract pressure 30 VAS onset 
                    elseif C{1,2}(i) == int_30(counter_30)
                     
                        onsets(b,3).onset = [onsets(b,3).onset;C{1,3}(i)/TR] ;
                        onsets(b,3).duration = [onsets(b,3).duration;((C{1,3}(i+1)/TR)-(C{1,3}(i)/TR))]; %subtracting stop  from start
                        counter_30 = counter_30 +1;
                        if counter_30 == 7
                            counter_30 = 6;
                        end
                end
            end

        end


    end
    end

end

%mkdir(['/projects/crunchie/nold/PEEP/fMRI/Data/' studyPart '/derivatives/spm_firstlevel/' sprintf('sub-%02.2d',subjects) filesep sprintf('ses-%02.2d',session)]);
%save(['/projects/crunchie/nold/PEEP/fMRI/Data/' studyPart '/derivatives/spm_firstlevel/' sprintf('sub-%02.2d',subjects) filesep sprintf('ses-%02.2d',session) filesep 'onsets.mat'],'onsets')

