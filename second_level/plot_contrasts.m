TR = 1.8;
bins = 12;

% mean
heat30_nacl = contrast.contrast(1:bins);
heat50_nacl = contrast.contrast(bins+1:2*bins);
heat70_nacl = contrast.contrast(2*bins+1:3*bins);

heat_nacl = contrast.contrast(1:3*bins);

pressure30_nacl = contrast.contrast(3*bins+1:4*bins);
pressure50_nacl = contrast.contrast(4*bins+1:5*bins);
pressure70_nacl = contrast.contrast(5*bins+1:6*bins);

pressure_nacl = contrast.contrast(3*bins+1:6*bins);

heat30_nax = contrast.contrast(6*bins+1:7*bins);
heat50_nax = contrast.contrast(7*bins+1:8*bins);
heat70_nax = contrast.contrast(8*bins+1:9*bins);

heat_nax = contrast.contrast(6*bins+1:9*bins);

pressure30_nax = contrast.contrast(9*bins+1:10*bins);
pressure50_nax = contrast.contrast(10*bins+1:11*bins);
pressure70_nax = contrast.contrast(11*bins+1:12*bins);

pressure_nax = contrast.contrast(9*bins+1:12*bins);

% se
heat30_nacl_se = contrast.standarderror(1:bins);
heat50_nacl_se = contrast.standarderror(bins+1:2*bins);
heat70_nacl_se = contrast.standarderror(2*bins+1:3*bins);

heat_nacl_se = contrast.standarderror(1:3*bins);

pressure30_nacl_se = contrast.standarderror(3*bins+1:4*bins);
pressure50_nacl_se = contrast.standarderror(4*bins+1:5*bins);
pressure70_nacl_se = contrast.standarderror(5*bins+1:6*bins);

pressure_nacl_se = contrast.standarderror(3*bins+1:6*bins);

heat30_nax_se = contrast.standarderror(6*bins+1:7*bins);
heat50_nax_se = contrast.standarderror(7*bins+1:8*bins);
heat70_nax_se = contrast.standarderror(8*bins+1:9*bins);

heat_nax_se = contrast.standarderror(6*bins+1:9*bins);

pressure30_nax_se = contrast.standarderror(9*bins+1:10*bins);
pressure50_nax_se = contrast.standarderror(10*bins+1:11*bins);
pressure70_nax_se = contrast.standarderror(11*bins+1:12*bins);

pressure_nax_se = contrast.standarderror(9*bins+1:12*bins);

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

%Nacl
% hFig = figure;
% subplot(2,1,1)
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p1 = plot((0:bins-1)*TR,heat30_nacl,'Color',grey,'LineWidth',2.5);hold on;
% p2 = plot((0:bins-1)*TR,heat50_nacl,'Color',blue,'LineWidth',2.5);
% p3 = plot((0:bins-1)*TR,heat70_nacl,'Color',red,'LineWidth',2.5);

% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:-1)*TR], [(heat30_nacl(2:bins-4)-heat30_nacl_se(2:bins-4)); flipud(heat30_nacl(2:bins-4)+heat30_nacl_se(2:bins-4))]',grey,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:-1)*TR], [(off(2:bins-4)-off_se(2:bins-4)); flipud(off(2:bins-4)+off_se(2:bins-4))]',blue,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% fi_p3 = fill([(0:bins-1)*TR (bins-1:-1:-1)*TR], [(off(2:bins-4)-off_se(2:bins-4)); flipud(off(2:bins-4)+off_se(2:bins-4))]',red,'EdgeColor', 'none');
% set(fi_p3, 'facealpha',0.3);
% 
% legend([p1 p2 p3],'30','50','70');
% xlabel('Time (s)','FontSize',12);
% ylabel('fMRI signal change (au)','FontSize',12);
% title('Heat Nacl');
% set(gca,'FontSize',8);
% box off
% 
% subplot(2,1,2)
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p4 = plot((0:bins-1)*TR,pressure30_nacl,'Color',grey,'LineWidth',2.5);hold on;
% p5 = plot((0:bins-1)*TR,pressure50_nacl,'Color',blue,'LineWidth',2.5);
% p6 = plot((0:bins-1)*TR,pressure70_nacl,'Color',red,'LineWidth',2.5);
% 
% xlabel('Time (s)','FontSize',12);
% ylabel('fMRI signal change (au)','FontSize',12);
% title('Pressure Nacl');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',8);
% box off
% % 
% % Naloxone
% hFig = figure;
% subplot(2,1,1)
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p7 = plot((0:bins-1)*TR,heat30_nax,'Color',grey,'LineWidth',2.5);hold on;
% p8 = plot((0:bins-1)*TR,heat50_nax,'Color',blue,'LineWidth',2.5);
% p9 = plot((0:bins-1)*TR,heat70_nax,'Color',red,'LineWidth',2.5);
% 
% xlabel('Time (s)','FontSize',12);
% ylabel('fMRI signal change (au)','FontSize',12);
% title('Heat Naloxone');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',8);
% box off
% 
% subplot(2,1,2)
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p10 = plot((0:bins-1)*TR,pressure30_nax,'Color',grey,'LineWidth',2.5);hold on;
% p11 = plot((0:bins-1)*TR,pressure50_nax,'Color',blue,'LineWidth',2.5);
% p12 = plot((0:bins-1)*TR,pressure70_nax,'Color',red,'LineWidth',2.5);
% 
% xlabel('Time (s)','FontSize',12);
% ylabel('fMRI signal change (au)','FontSize',12);
% title('Pressure Naloxone');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',8);
% box off


% Heat Nacl vs. Heat Naloxone (50 VAS)
% hFig = figure;
% subplot(2,1,1);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p13 = plot((0:bins-1)*TR,heat50_nacl,'Color',blue2,'LineWidth',2.5);hold on;
% p14 = plot((0:bins-1)*TR,heat50_nax,'Color',orange,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% xline(1,'Color',grey,'LineWidth',2);
% xline(5,'Color',red,'LineWidth',2,'LineStyle','-');
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat50_nacl(1:bins)-heat50_nacl_se(1:bins)); flipud(heat50_nacl(1:bins)+heat50_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat50_nax(1:bins)-heat50_nax_se(1:bins)); flipud(heat50_nax(1:bins)+heat50_nax_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('Heat [VAS 50]');
% legend([p13 p14],'SAL','NLX');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% Pressure Nacl vs. Heat Naloxone (70 VAS)
% hfig = figure;
% subplot(2,1,2);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p15 = plot((0:bins-1)*TR,pressure50_nacl,'Color',blue2,'LineWidth',2.5);hold on;
% p16 = plot((0:bins-1)*TR,pressure50_nax,'Color',orange,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% xline(1,'Color',grey,'LineWidth',2);
% xline(3,'Color',red,'LineWidth',2,'LineStyle','-');
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure50_nacl(1:bins)-pressure50_nacl_se(1:bins)); flipud(pressure50_nacl(1:bins)+pressure50_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure50_nax(1:bins)-pressure50_nax_se(1:bins)); flipud(pressure50_nax(1:bins)+pressure50_nax_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('Pressure [VAS 70]');
% legend([p15 p16],'SAL','NLX');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% % xlim([-1 46]);
% set(gcf,'PaperPositionMode','auto');
% time = (0:bins-1)*TR;
% time = time(:);
%tbl = table(time,heat70_nacl,heat70_nax,pressure70_nacl,pressure70_nax,heat70_nacl_se,heat70_nax_se,pressure70_nacl_se,pressure70_nax_se);
%writetable(tbl,'/home/nold/Desktop/PEEP/fMRI/Figures/FIR_model_treat_gender_heat_OPERC.csv');

% Heat Nacl vs. Heat Naloxone (70 VAS)
% hFig = figure;
% %subplot(2,1,1);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p13 = plot((0:bins-1)*TR,heat70_nacl,'Color',blue2,'LineWidth',2.5);hold on;
% p14 = plot((0:bins-1)*TR,heat70_nax,'Color',orange,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% xline(1,'Color',grey,'LineWidth',2);
% xline(5,'Color',red,'LineWidth',2,'LineStyle','-');
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nacl(1:bins)-heat70_nacl_se(1:bins)); flipud(heat70_nacl(1:bins)+heat70_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nax(1:bins)-heat70_nax_se(1:bins)); flipud(heat70_nax(1:bins)+heat70_nax_se(1:bins))]',orange,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('Heat [VAS 70]');
% legend([p13 p14],'SAL','NLX');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% % xlim([-1 46]);
% set(gcf,'PaperPositionMode','auto');
% time = (0:bins-1)*TR;
% time = time(:);
% 
% tbl = table(time,heat70_nacl,heat70_nax,heat70_nacl_se,heat70_nax_se);
% writetable(tbl,'/home/nold/Desktop/PEEP/fMRI_new/betas/FIR_inter_treat_int_heat_pag_-2_-24_-8.csv');

% %% Pressure Nacl vs. Heat Saline (70 VAS)
% hFig = figure;
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p15 = plot((0:bins-1)*TR,pressure70_nacl,'Color',blue2,'LineWidth',2.5);hold on;
% p16 = plot((0:bins-1)*TR,heat70_nacl,'Color',orange,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% xline(1,'Color',grey,'LineWidth',2);
% xline(3,'Color',red,'LineWidth',2,'LineStyle','-');
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nacl(1:bins)-pressure70_nacl_se(1:bins)); flipud(pressure70_nacl(1:bins)+pressure70_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nacl(1:bins)-heat70_nacl_se(1:bins)); flipud(heat70_nacl(1:bins)+heat70_nacl_se(1:bins))]',orange,'EdgeColor', 'none');
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
%tbl = table(time,heat70_nacl,heat70_nax,pressure70_nacl,pressure70_nax,heat70_nacl_se,heat70_nax_se,pressure70_nacl_se,pressure70_nax_se);
%writetable(tbl,'/home/nold/Desktop/PEEP/fMRI/Figures/FIR_model_treat_gender_heat_OPERC.csv');



% %% Heat vs Pressure NLX
% 
% % Heat Nacl vs. Heat Naloxone (70 VAS)
% hFig = figure;
% subplot(2,1,1);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p13 = plot((0:bins-1)*TR,pressure70_nax,'Color',grey2,'LineWidth',2.5);hold on;
% p14 = plot((0:bins-1)*TR,heat70_nax,'Color',black,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nax(1:bins)-pressure70_nax_se(1:bins)); flipud(pressure70_nax(1:bins)+pressure70_nax_se(1:bins))]',grey2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nax(1:bins)-heat70_nax_se(1:bins)); flipud(heat70_nax(1:bins)+heat70_nax_se(1:bins))]',black,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('NLX');
% legend([p13 p14],'Pressure70','Heat70');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% % Pressure Nacl vs. Heat Nacl (70 VAS)
% subplot(2,1,2);
% set(hFig,'units','pixel','pos',[900 400 500 350]);
% p15 = plot((0:bins-1)*TR,pressure70_nacl,'Color',grey2,'LineWidth',2.5);hold on;
% p16 = plot((0:bins-1)*TR,heat70_nacl,'Color',black,'LineWidth',2.5);
% yline(0,'Color',grey,'LineWidth',2);
% 
% fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nacl(1:bins)-pressure70_nacl_se(1:bins)); flipud(pressure70_nacl(1:bins)+pressure70_nacl_se(1:bins))]',grey2,'EdgeColor', 'none');
% set(fi_p1, 'facealpha',0.3);
% fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(heat70_nacl(1:bins)-heat70_nacl_se(1:bins)); flipud(heat70_nacl(1:bins)+heat70_nacl_se(1:bins))]',black,'EdgeColor', 'none');
% set(fi_p2, 'facealpha',0.3);
% 
% xlabel('Time (s)','FontSize',20);
% ylabel('fMRI signal change (au)','FontSize',20);
% title('SAL');
% legend([p15 p16],'Pressure70','Heat70');
% 
% set(gcf,'color','w');
% set(gca,'FontSize',15);
% box off
% 
% % xlim([-1 46]);
% set(gcf,'PaperPositionMode','auto');


%% Pressure Nacl 70 vs. Pressure NLX (70 VAS)
hFig = figure;
set(hFig,'units','pixel','pos',[900 400 500 350]);
p15 = plot((0:bins-1)*TR,pressure70_nacl,'Color',blue2,'LineWidth',2.5);hold on;
p16 = plot((0:bins-1)*TR,pressure70_nax,'Color',orange,'LineWidth',2.5);
yline(0,'Color',grey,'LineWidth',2);
xline(1,'Color',grey,'LineWidth',2);
xline(3,'Color',red,'LineWidth',2,'LineStyle','-');

fi_p1 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nacl(1:bins)-pressure70_nacl_se(1:bins)); flipud(pressure70_nacl(1:bins)+pressure70_nacl_se(1:bins))]',blue2,'EdgeColor', 'none');
set(fi_p1, 'facealpha',0.3);
fi_p2 = fill([(0:bins-1)*TR (bins-1:-1:0)*TR], [(pressure70_nax(1:bins)-pressure70_nax_se(1:bins)); flipud(pressure70_nax(1:bins)+pressure70_nax_se(1:bins))]',orange,'EdgeColor', 'none');
set(fi_p2, 'facealpha',0.3);

xlabel('Time (s)','FontSize',20);
ylabel('fMRI signal change (au)','FontSize',20);
title('');
legend([p15 p16],'Saline','Naloxone');

set(gcf,'color','w');
set(gca,'FontSize',15);
box off

% xlim([-1 46]);
set(gcf,'PaperPositionMode','auto');
time = (0:bins-1)*TR;
time = time(:);
tbl = table(time,heat70_nacl,heat70_nax,pressure70_nacl,pressure70_nax,heat70_nacl_se,heat70_nax_se,pressure70_nacl_se,pressure70_nax_se);
writetable(tbl,'/home/nold/Desktop/PEEP/fMRI_new/betas/heat_pressure/FIR_model_treat_int_pressure_36_3_18.csv');
