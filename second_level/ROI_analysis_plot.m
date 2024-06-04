
nlx = heat30_exHi_nlx
sal = heat30_exHi_sal
mod_id = 2;
int = 30;
ex_int = 1;


% load covariate (FTP value z-standardised)
fid = fopen ('/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/participants.tsv');
C = textscan (fid, '%s %d %d %d %d %d %f %f', 'HeaderLines', 1);
fclose (fid);
pwc_vec = C{8};
gender_vec = double(C{3});

xfs = 10;
ms = 12;
lw = 2;
col.nacl = [2 72 115]/255;
col.nlx = [242 149 68]/255;

[h,p,ci,stats]=ttest(sal,nlx)

% Correlation
[r,p,rl,ru] = corrcoef(pwc_vec,nlx)
[r,p,rl,ru] = corrcoef(pwc_vec,sal)


% hFig = figure;
% set(hFig,'units','pixel','pos',[900 400 [185*2 135*2]]);
% plot(pwc_vec,sal,'.','MarkerSize',ms,'MarkerFaceColor',col.nacl,'MarkerEdgeColor',col.nacl);hold on;
% plot(pwc_vec,nlx,'.','MarkerSize',ms,'MarkerFaceColor',col.nlx,'MarkerEdgeColor',col.nlx);
% h = lsline;
% set(h(2),'LineWidth',lw,'color',col.nacl);
% set(h(1),'LineWidth',lw,'color',col.nlx);
% legend('sal','nlx')
% xlim([0.3 3.5]);
% %ylim([-3 3]);
% ylabel('fMRI signal (au)','FontSize',10);
% xlabel('FTP (Watt/kg)','FontSize',10);
% box off
% set(gcf,'color','w');
% set(gcf,'PaperPositionMode','auto');
% title('Overall')


% Calculate Linear Rwegression with Interaction Effects 
betas = [sal;nlx];
pharm_cond = [-ones(length(sal),1);ones(length(nlx),1)];
modality = [repelem(mod_id,length(sal))';repelem(mod_id,length(sal))'];
gender = [gender_vec;gender_vec];
int_vec = repelem(int,length(nlx)*2)';
ex_int_vec = repelem(ex_int,length(nlx)*2)';
ftp = [pwc_vec;pwc_vec];
%VAS = repelem(int,length(ftp),1);

tbl = table(pharm_cond,ftp,gender,modality,int_vec,ex_int_vec,betas,'VariableNames',{'pharm','ftp','gender','modality','int','ex_int','betas'});
%mdl = fitlm(tbl,'betas ~ ftp + gender + pharm + gender*ftp*pharm')
%mdl = fitlm(tbl,'betas ~ modality + pharm + modality*pharm')
writetable(tbl,'/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_secondlevel/roi/roi_betas_heat30_exHi_frontMidline_vega.csv');




%% ----------- Male only
% sal_male= sal(gender_vec==-1);
% nlx_male = nlx(gender_vec==-1);
% 
% % T-test
% [h,p,ci,stats]=ttest(sal_male,nlx_male)
% 
% % Correlation
% [r,p,rl,ru] = corrcoef(pwc_vec(gender_vec ==-1),nlx_male)
% [r,p,rl,ru] = corrcoef(pwc_vec(gender_vec ==-1),sal_male)
% 
% %% Calculate Linear Regression with Interaction Effects 
% betas = [sal_male;nlx_male];
% pharm_cond = [-ones(length(sal_male),1);ones(length(nlx_male),1)];
% ftp = [pwc_vec(gender_vec==-1);pwc_vec(gender_vec==-1)];
% 
% tbl = table(ftp,pharm_cond,betas,'VariableNames',{'ftp','pharm','betas'});
% mdl = fitlm(tbl,'betas ~ ftp*pharm')
% 
% hFig = figure;
% set(hFig,'units','pixel','pos',[900 400 [185*2 135*2]]);
% plot(pwc_vec(gender_vec==-1),sal_male,'.','MarkerSize',ms,'MarkerFaceColor',col.nacl,'MarkerEdgeColor',col.nacl);hold on;
% plot(pwc_vec(gender_vec==-1),nlx_male,'.','MarkerSize',ms,'MarkerFaceColor',col.nlx,'MarkerEdgeColor',col.nlx);
% h = lsline;
% set(h(2),'LineWidth',lw,'color',col.nacl);
% set(h(1),'LineWidth',lw,'color',col.nlx);
% legend('sal','nlx')
% xlim([0.3 3.5]);
% %ylim([-3 3]);
% ylabel('fMRI signal (au)','FontSize',xfs);
% xlabel('FTP (Watt/kg)','FontSize',xfs);
% box off
% set(gcf,'color','w');
% set(gcf,'PaperPositionMode','auto');
% title('Male');
% 
% %% ----------- Female only
% sal_female = sal(gender_vec==1);
% nlx_female = nlx(gender_vec==1);
% 
% % T test
% [h,p,ci,stats]=ttest(sal_female,nlx_female)
% 
% % Correlation
% [r,p,rl,ru] = corrcoef(pwc_vec(gender_vec ==1),nlx_female)
% [r,p,rl,ru] = corrcoef(pwc_vec(gender_vec ==1),sal_female)
% 
% %% Calculate Linear Regression with Interaction Effects 
% betas = [sal_female;nlx_female];
% pharm_cond = [-ones(length(sal_female),1);ones(length(nlx_female),1)];
% ftp = [pwc_vec(gender_vec==1);pwc_vec(gender_vec==1)];
% 
% tbl = table(ftp,pharm_cond,betas,'VariableNames',{'ftp','pharm','betas'});
% mdl = fitlm(tbl,'betas ~ ftp + pharm + ftp*pharm')
% 
% hFig = figure;
% set(hFig,'units','pixel','pos',[900 400 [185*2 135*2]]);
% plot(pwc_vec(gender_vec==1),sal_female,'.','MarkerSize',ms,'MarkerFaceColor',col.nacl,'MarkerEdgeColor',col.nacl);hold on;
% plot(pwc_vec(gender_vec==1),nlx_female,'.','MarkerSize',ms,'MarkerFaceColor',col.nlx,'MarkerEdgeColor',col.nlx);
% h = lsline;
% set(h(2),'LineWidth',lw,'color',col.nacl);
% set(h(1),'LineWidth',lw,'color',col.nlx);
% legend('sal','nlx')
% xlim([0.3 3.5]);
% %ylim([-3 3]);
% ylabel('fMRI signal (au)','FontSize',xfs);
% xlabel('FTP (Watt/kg)','FontSize',xfs);
% box off
% set(gcf,'color','w');
% set(gcf,'PaperPositionMode','auto');
% title('Female')
% 
% 
% %% T test differences gender
% % T test
% [h,p,ci,stats]=ttest2(sal_male,sal_female)
% % T test
% [h,p,ci,stats]=ttest2(nlx_male,nlx_female)
% 
% 
% %% Barplot
% betas = [sal,nlx];
% ms = 10;
% size1plot = [185*2 135*2];
% afs = 10;
% xfs = 15;
% yfs = 15;
% bw = 0.8;
% dev = 0.14;
% err = 1; % ci 90%: 1.645
% sr1 = sqrt(length(gender_vec(gender_vec==1)));
% sr2 = sqrt(length(gender_vec(gender_vec==-1)));
% 
% val_bar = [mean(betas(gender_vec==1,1)) mean(betas(gender_vec==1,2));...
%     mean(betas(gender_vec==-1,1)) mean(betas(gender_vec==-1,2))];
% val_se = [std(betas(gender_vec==1,1))/sr1 * err std(betas(gender_vec==1,2))/sr1 * err;...
%     std(betas(gender_vec==-1,1))/sr2 * err std(betas(gender_vec==-1,2))/sr2 * err];
% 
% hFig = figure;
% set(hFig,'units','pixel','pos',[900 400 size1plot]);
% h = bar(val_bar,'EdgeColor','none','BarWidth',bw); hold on;
% set(h(1),'FaceColor',col.nacl);
% set(h(2),'FaceColor',col.nlx); xpos = [1-dev 1+dev 2-dev 2+dev];
% y = [val_bar(1,:)-val_se(1,:) val_bar(2,:)-val_se(2,:); val_bar(1,:)+val_se(1,:) val_bar(2,:)+val_se(2,:)];
% plot([xpos; xpos], y, '-k', 'LineWidth',lw);
% plot([xpos; xpos], [val_bar(1,:) val_bar(2,:)], '.k', 'MarkerSize',ms-4);
% % errorbar([1-dev 1+dev; 2-dev 2+dev; 3-dev 3+dev],val_bar,val_se,'.','Color','k', 'LineWidth',lw);
% set(gca,'xtick',[1 2],'xticklabel',{'Female','Male'});
% ax = gca;
% ax.XAxis.FontSize = xfs;
% ax.YAxis.FontSize = afs;
% ylabel('fMRI signal change (au)','FontSize',yfs);
% box off
% %legend('SAL','NLX');
% title(['Coordinate: ' num2str(coordinate)]);
% set(gcf,'color','w');
% set(gcf,'PaperPositionMode','auto');