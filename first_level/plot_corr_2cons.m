clear
subs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 
exclude    = [8];

if ~isempty(exclude)
    subs = subs(~ismember(subs,exclude));
end

base_dir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_firstlevel/HRF_complete_onset_stimulus_6mm/';
coordinate = [39 -8 12];

% Heat SAL 
% neg Cor FTP
% 62 -30 -15
% -69 -23 -9
% 6 45 9
% 42 -27 15

% Heat NLX:
% neg Cor FTP
% 2 -14 5
% -16 63 -4
% -8 32 -8

%pos cor
% 39 -24 36




pub = 1;
image_nr = [27 31]; 

% heat sal 3
% heat NLX 4
% Pressure sal 5
% pressure nlx 6
%Heat NLX 21 22 % Heat SAL 19 23 pressure SAL 23 24 %pressure NLX 25 26

betas = zeros(numel(subs),numel(image_nr));

for g = 1:size(subs,2)
    name = sprintf('sub-%02.2d',subs(g));
    subdir = fullfile(base_dir,name);

    cons = spm_select('FPList',subdir,'^s6w_nlco_dartelcon.*.nii$');
    con = cons(image_nr,:);

    for j = 1:size(con,1)
        V = spm_vol(con(j,:));

        if g == 1 && j == 1
            v2m = spm_get_space(V(1).fname);
            m2v = inv(v2m);

            for i=1:size(coordinate,1)
                vox(i,1:3) = coordinate(i,:)*m2v(1:3,1:3) + m2v(1:3,4)';
            end

            disp(vox)

        end
        betas(g,j) = spm_get_data(V,vox');

    end

end

fid = fopen ('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/participants.tsv');
C = textscan (fid, '%s %d %s %d %d %d %f %f', 'HeaderLines', 1);
fclose (fid);
ftp_vec = double(C{5});
pwc_vec = double(C{8});
%gender_vec = double(C{3});


%betas = betas(gender_vec==1,:);

col.low = [27 90 129]/255;
col.high = [243 159 86]/255;

if pub
    lw = 2;
    ms = 20;
    sizecorr = [185*2 135*2]; %for figures and publication
    afs = 12;
    xfs = 15;
    yfs = 15;
else
    lw = 1.5;
    ms = 15;
    sizecorr = [320 250];
    afs = 10;
    xfs = 12;
    yfs = 12;
end

xlims = [0 3.5];
ylims = [5*floor(min(betas-1)/5) 5*ceil(max(betas+1)/5)];

hFig = figure;
set(hFig,'units','pixel','pos',[900 400 sizecorr]);

plot(pwc_vec,betas,'.');hold on;

p1 = plot(pwc_vec,betas(:,1),'.','MarkerFaceColor', col.low,'MarkerEdgeColor',col.low,'MarkerSize',ms);hold on
h = lsline;
set(h(1),'color',col.low,'LineWidth',lw);

p2 = plot(pwc_vec,betas(:,2),'.','MarkerFaceColor', col.high,'MarkerEdgeColor',col.high,'MarkerSize',ms); hold on
% p3 = plot(ratings(r.groups==3),-betas(r.groups==3),'.','MarkerFaceColor', col.r100,'MarkerEdgeColor',col.r100,'MarkerSize',ms);
h = lsline;
set(h(1),'color', col.high,'LineWidth',lw);

set(gca,'FontSize',afs);

xlabel('relative FTP (FTP/kg)','FontSize',xfs);
ylabel(['fMRI signal change (\Deltaau)'],'FontSize',yfs);
xlim(xlims);
%ylim(ylims);

box off
set(gcf,'color','w');

legend('Exercise Low','Exercise High');

set(gcf,'PaperPositionMode','auto');
title(['Coordinate: ' num2str(coordinate)]);


cd('/home/nold/Desktop/PEEP/fMRI_new/betas/');
%exportgraphics(gcf,sprintf('cor_2cons_exercise_diff_%d_%d_%d.png',coordinate(1),coordinate(2),coordinate(3)),'Resolution',500);

tbl = table(betas);
writetable(tbl,sprintf('cor_betas_exercise_heat_nlx_%d_%d_%d.csv',coordinate(1),coordinate(2),coordinate(3)));

% rh = corr(ratings, -betas);
% 
% if rh < 0
%     tpos = [xlims(1) ylims(1)]/1.2;
% %     legend([p1 p2 p3],'NaCl','Remi50','Remi100','location','northeast');
% elseif rh > 0
%     tpos = [xlims(2)-20 ylims(1)]/1.2;
% %     legend([p1 p2 p3],'NaCl','Remi50','Remi100','location','northwest');
% end
% text(tpos(1),tpos(2),['r = ' num2str(sprintf('%.2f',rh))],'FontSize',afs,'color',[0 0 0],'FontWeight','bold');  %


% print -dsvg -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_57_6_3.svg
%print -dsvg -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_6_51_-6.svg
% print -dtiff -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_6_6_3.tiff
% print -dtiff -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_8_50_-6.tiff




