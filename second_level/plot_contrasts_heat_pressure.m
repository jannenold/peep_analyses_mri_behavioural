TR = 1.8;
bins = 12;

% mean
heat30_nacl_lo = contrast.contrast(1:bins);
heat50_nacl_lo = contrast.contrast(bins+1:2*bins);
heat70_nacl_lo = contrast.contrast(2*bins+1:3*bins);

heat_nacl_lo = mean([heat30_nacl_lo,heat50_nacl_lo,heat70_nacl_lo],2);

pressure30_nacl_lo = contrast.contrast(3*bins+1:4*bins);
pressure50_nacl_lo = contrast.contrast(4*bins+1:5*bins);
pressure70_nacl_lo = contrast.contrast(5*bins+1:6*bins);

pressure_nacl_lo = mean([pressure30_nacl_lo,pressure50_nacl_lo,pressure70_nacl_lo],2);

% se
heat30_nacl_se = contrast.standarderror(1:bins);
heat50_nacl_se = contrast.standarderror(bins+1:2*bins);
heat70_nacl_se = contrast.standarderror(2*bins+1:3*bins);

heat_nacl_se = contrast.standarderror(1:3*bins);

pressure30_nacl_se = contrast.standarderror(3*bins+1:4*bins);
pressure50_nacl_se = contrast.standarderror(4*bins+1:5*bins);
pressure70_nacl_se = contrast.standarderror(5*bins+1:6*bins);

pressure_nacl_se = contrast.standarderror(3*bins+1:6*bins);

% colours 
red = [205/255 0/255 0/255];
blue = [24/255 116/255 205/255];
grey = [0.7 0.7 0.7];
grey2 = [0.5 0.5 0.5];
black = [1/255 1/255 1/255];
col  = [5 4 4];

blue2 = [2 72 115]/255;
orange = [242 149 68]/255;

% =========
% Plots
%==========

% %% Pressure Nacl vs. Heat Saline 
% hFig = figure;
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p15 = plot((0:bins-1)*TR,pressure_nacl_lo,'Color',blue2,'LineWidth',2.5);hold on;
% p16 = plot((0:bins-1)*TR,heat_nacl_lo,'Color',orange,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure_nacl_lo(1:bins)-pressure_nacl_se(1:bins)); flipud(pressure_nacl_lo(1:bins)+pressure_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat_nacl_lo(1:bins)-heat_nacl_se(1:bins)); flipud(heat_nacl_lo(1:bins)+heat_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('Saline');
% legend([p15 p16],'Pressure','Heat');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% % xlim([-1 46]);
% set(gcf,'PaperPositionMode','auto');
% time = (0:bins-1)*TR;
% time = time(:);
% tbl = table(time,pressure_nacl_lo,heat_nacl_lo,pressure_nacl_se(1:12),heat_nacl_se(1:12));
% writetable(tbl,'/home/nold/Desktop/PEEP/fMRI_new/betas/heat_pressure/FIR/FIR_model_heat_pressure_-22_-62_57.csv');


%% Pressure Nacl vs. Heat Saline (70 VAS)
hFig = figure;
set(hFig,'units','pixel','pos',[900 400 500 350]);
p15 = plot((0:bins-1)*TR,pressure70_nacl_lo,'Color',blue2,'LineWidth',2.5);hold on;
p16 = plot((0:bins-1)*TR,heat70_nacl_lo,'Color',orange,'LineWidth',2.5);
yline(0,'Color',grey,'LineWidth',2);
xline(1,'Color',grey,'LineWidth',2);
xline(3,'Color',red,'LineWidth',2,'LineStyle','-');

fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nacl_lo(1:bins)-pressure70_nacl_se(1:bins)); flipud(pressure70_nacl_lo(1:bins)+pressure70_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
set(fi_p1, 'facealpha',0.3);
fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nacl_lo(1:bins)-heat70_nacl_se(1:bins)); flipud(heat70_nacl_lo(1:bins)+heat70_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
set(fi_p2, 'facealpha',0.3);

xlabel('Time (s)','FontSize',20);
ylabel('fMRI signal change (au)','FontSize',20);
title('Saline');
legend([p15 p16],'Pressure','Heat');

set(gcf,'color','w');
set(gca,'FontSize',15);
box off

% xlim([-1 46]);
set(gcf,'PaperPositionMode','auto');
time = (0:bins-1)*TR;
time = time(:);
% %tbl = table(time,heat70_nacl,heat70_nax,pressure70_nacl,pressure70_nax,heat70_nacl_se,heat70_nax_se,pressure70_nacl_se,pressure70_nax_se);
% %writetable(tbl,'/home/nold/Desktop/PEEP/fMRI/Figures/FIR_model_treat_gender_heat_OPERC.csv');
% 
%% Parametric Effect Heat (30,50,70 VAS)
% 
% hFig = figure;
% subplot(2,1,1);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p15 = plot((0:bins-1)*TR,heat30_nacl_lo,'Color',blue2,'LineWidth',2.5);hold on;
% p16 = plot((0:bins-1)*TR,heat50_nacl_lo,'Color',orange,'LineWidth',2.5);
% p17 = plot((0:bins-1)*TR,heat70_nacl_lo,'Color',grey,'LineWidth',2.5);
% 
% yline(0,'Color',grey,'LineWidth',2);
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat30_nacl_lo(1:bins)-heat30_nacl_se(1:bins)); flipud(heat30_nacl_lo(1:bins)+heat30_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat50_nacl_lo(1:bins)-heat50_nacl_se(1:bins)); flipud(heat50_nacl_lo(1:bins)+heat50_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nacl_lo(1:bins)-heat70_nacl_se(1:bins)); flipud(heat70_nacl_lo(1:bins)+heat70_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('Heat');
% legend([p15 p16 p17],'30','50','70');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% set(gcf,'PaperPositionMode','auto');
% time = (0:bins-1)*TR;
% time = time(:);
% 
% subplot(2,1,2);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p15 = plot((0:bins-1)*TR,pressure30_nacl_lo,'Color',blue2,'LineWidth',2.5);hold on;
% p16 = plot((0:bins-1)*TR,pressure50_nacl_lo,'Color',orange,'LineWidth',2.5);
% p17 = plot((0:bins-1)*TR,pressure70_nacl_lo,'Color',grey,'LineWidth',2.5);
% 
% yline(0,'Color',grey,'LineWidth',2);
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure30_nacl_lo(1:bins)-pressure30_nacl_se(1:bins)); flipud(pressure30_nacl_lo(1:bins)+pressure30_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure50_nacl_lo(1:bins)-pressure50_nacl_se(1:bins)); flipud(pressure50_nacl_lo(1:bins)+pressure50_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nacl_lo(1:bins)-pressure70_nacl_se(1:bins)); flipud(pressure70_nacl_lo(1:bins)+pressure70_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('Pressure');
% legend([p15 p16 p17],'30','50','70');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% set(gcf,'PaperPositionMode','auto');
% time = (0:bins-1)*TR;
% time = time(:);

% tbl = table(time,heat30_nacl_lo,heat50_nacl_lo,heat70_nacl_lo,heat30_nacl_se,heat50_nacl_se,heat70_nacl_se,...
%    pressure30_nacl_lo,pressure50_nacl_lo,pressure70_nacl_lo,pressure30_nacl_se,pressure50_nacl_se,pressure70_nacl_se);
% writetable(tbl,'/home/nold/Desktop/PEEP/fMRI_new/betas/heat_pressure/FIR/FIR_model_heat_pressure_inter_mod_int_-28_-16_56.csv');
% 

