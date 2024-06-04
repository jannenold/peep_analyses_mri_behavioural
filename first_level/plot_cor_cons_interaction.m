clear
subs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 
exclude    = [8];

if ~isempty(exclude)
    subs = subs(~ismember(subs,exclude));
end

base_dir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_firstlevel/HRF_complete_onset_stimulus_6mm/';


coordinate = {[-14 62 -4]};


% Saline heat:
% 6 45 10

% Naloxone Heaz
% -14 62 -4


image_nr = [27 31];

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
C = textscan (fid, '%s %d %d %d %d %d %f %f', 'HeaderLines', 1);
fclose (fid);
ftp_vec = double(C{5});
pwc_vec = double(C{8});
gender_vec = double(C{3});

betas = betas(gender_vec==1,:);
pwc_vec = pwc_vec(gender_vec==1,:);
pwc_vec_female = pwc_vec(gender_vec==1);
pwc_vec_male = pwc_vec(gender_vec==-1);
pwc_vec_new = pwc_vec_female;


col.nacl = [2 72 115]/255;
col.nlx = [242 149 68]/255;
col.re50 = [68 183 194]/255;
col.r100 = [2 75 122]/255;


if pub
    lw = 2;
    ms = 25;
    sizecorr = [185*2 135*2]; %for figures and publication
    afs = 12;
    xfs = 15;
    yfs = 15;
else
    lw = 1.5;
    ms = 25;
    sizecorr = [320 250];
    afs = 10;
    xfs = 12;
    yfs = 12;
end

xlims = [0 4];
ylims = [5*floor(min(betas-1)/5) 5*ceil(max(betas+1)/5)];

hFig = figure;
set(hFig,'units','pixel','pos',[900 400 sizecorr]);

plot(pwc_vec_new,betas,'.');hold on;
h = lsline;
set(h(1),'color',[0 0 0],'LineWidth',lw);

p1 = plot(pwc_vec_female,betas,'.','MarkerFaceColor', col.nacl,'MarkerEdgeColor',col.nacl,'MarkerSize',ms);hold on
p2 = plot(pwc_vec,betas,'.','MarkerFaceColor', col.nlx,'MarkerEdgeColor',col.nlx,'MarkerSize',ms);hold on
% p3 = plot(ratings(r.groups==3),-betas(r.groups==3),'.','MarkerFaceColor', col.r100,'MarkerEdgeColor',col.r100,'MarkerSize',ms);

set(gca,'FontSize',afs);

xlabel('relative FTP [Watt/kg]','FontSize',xfs);
ylabel(['fMRI signal change (\Deltaau)'],'FontSize',yfs);
xlim(xlims);
ylim(ylims);

box off
set(gcf,'color','w');

set(gcf,'PaperPositionMode','auto');

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
coordinate = coordinate{:};
%cd('/home/nold/Desktop/PEEP/fMRI/Figures/');
%exportgraphics(gcf,sprintf('cor_cons_exercise_%d_%d_%d.png',coordinate(1),coordinate(2),coordinate(3)),'Resolution',500);


% print -dsvg -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_57_6_3.svg
%print -dsvg -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_6_51_-6.svg
% print -dtiff -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_6_6_3.tiff
% print -dtiff -painters -noui -r500 /home/tinnermann/remi3/Figures/Corr_rating_8_50_-6.tiff




