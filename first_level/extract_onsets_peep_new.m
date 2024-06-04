function[onsets] = extract_onsets_peep_new(subjects,session,studyPart,nr,model)

dataDir = '/projects/crunchie/nold/PEEP/behavioural';
blocks = nr;
TR = 1.8; % in seconds

for b = 1:numel(blocks)
	
		C = [];
		D = [];
		
        % get intensities
		fid = fopen([dataDir,filesep,studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_behav.tsv',subjects,blocks(b),session+1)]);
		D = textscan(fid,'%f %f %f %f %s %f %f %f %f %f %f','HeaderLines',1);
		fclose(fid);
		
        % get onsets
		fid = fopen([dataDir,filesep,studyPart,filesep,sprintf('sub-%02.2d/ses-%02.2d/',subjects,session),sprintf('sub-%02.2d_block-%02.2d-peep_exp-day-%02.2d_all_events.tsv',subjects,blocks(b),session+1)]);
		C = textscan(fid,'%f %f %f %s','HeaderLines',1);
		fclose(fid);
		
        % combine intensities & onsets
		tabC = table(repmat(b,size(C{:,3},1),1),C{:,3},C{:,3}./TR,C{:,4},'VariableNames',{'Block','t','Scan','Event'});
		relIx = cell2mat(cellfun(@(x) ~isempty(regexp(x,'^start_plateau','ONCE')),tabC.Event,'UniformOutput',false)); % get relevant indices
		tabC(~relIx,:) = [];
		tabC.Modality = double(strcmp(tabC.Event,'start_plateau_pressure'));
		tabC.nModality = categorical(tabC.Modality,[1 0],{'Pressure','Heat'});
		tabC.Intensity = D{1,6};
		
        % onset extraction
		modalities = {'pressure','heat'};
		intensities = [70 50 30];
		onsetsX = struct;
		xc = 0;

		for xm = 1:numel(modalities)
			for xii = 1:numel(intensities)
				xc = xc + 1;
				onsetsX(b,xc).name = sprintf('%s_%d',modalities{xm},intensities(xii));
				onsetsX(b,xc).onset = tabC.Scan(tabC.nModality==modalities{xm} & tabC.Intensity==intensities(xii));
				
                if strcmp(model,'FIR')
                onsetsX(b,xc).duration = 0;

                elseif strcmp(model,'HRF')
                    onsetsX(b,xc).duration = 0;
                    
                    onsetsX(b,xc).tmod = 0;
                    onsetsX(b,xc).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    onsetsX(b,xc).orth = 1;

                end 
			end
		end
		
        onsets = onsetsX;
end