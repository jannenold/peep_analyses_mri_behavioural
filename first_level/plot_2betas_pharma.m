function plot_2betas_pharma

subs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48]; 
exclude    = [8];

if ~isempty(exclude)
    subs = subs(~ismember(subs,exclude));
end

base_dir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/derivatives/spm_firstlevel/HRF_complete_onset_stimulus_6mm';

barplot = 1;
%plot individual data points
ind_dots = 0;
%figure size publication or not
pub = 1;

coordinate = [46 8 6];
image_nr = [3 4];


% Heat NLX > SAL
% [46 8 6];
% [3 -22 -9];


% Heat SAL > NLX
%[46 -4 0]; planum polare
% [-26 -16.5 -30]; %amygdala

% Pressure NLX > SAL
% [50 -21 62]; %postc. gyrus
% [9 -10 -3]; % thalamus
% [30 -16.5 -30]; % parahippo. gyrus

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
ftp_vec = double(C{5});
ftp_median = 122;
gender_vec = C{3};


col.base = [128 175 205]/255;
col.infs = [51 107 156]/255;

col.nacl = [45 154 166]/255;
col.nlx = [242 92 5]/255;

%col.nacl = [0.7 0.7 0.7];


bw = 0.8;
dev = 0.14;
err = 1; % ci 90%: 1.645

if pub
    lw = 1;
    size1plot = [165 135]; %for figures and publication
    afs = 6;
    xfs = 7.5;
    yfs = 8;
else
    lw = 1;
    ms = 15;
    size1plot = [320 250];
    afs = 10;
    xfs = 11;
    yfs = 12;
end

if barplot == 1
    ms = 10;
    size1plot = [200 135*3]; 
    afs = 10;
    xfs = 15;
    yfs = 15;
    sr1 = sqrt(length(ftp_vec));
    sr2 = sqrt(length(ftp_vec));
    
    val_bar = [mean(betas(:,1)) mean(betas(:,2))];
    val_se = [std(betas(:,1))/sr1 * err std(betas(:,2))/sr1 * err];
    
    hFig = figure;
    set(hFig,'units','pixel','pos',[900 400 size1plot]);
    
    h = bar(val_bar,'EdgeColor','none','BarWidth',bw,'FaceColor','flat'); hold on;

    h(1).CData(1,:) = col.nacl;
    h(1).CData(2,:) = col.nlx;

    %set(h,'FaceColor',col.nacl);
    %set(b,'FaceColor',col.nlx);
    
    if ind_dots
        set(h,'FaceAlpha',0.6)
        a = -0.05;
        b = 0.05;
        plot([1-dev+((b-a).*rand(1,numel(betas(:,1)))+a)],betas(:,1),'.','MarkerFaceColor', col.nacl,'MarkerEdgeColor',col.nacl,'MarkerSize',ms);
        plot([1+dev+((b-a).*rand(1,numel(betas(:,2)))+a)],betas(:,2),'.','MarkerFaceColor', col.nlx,'MarkerEdgeColor',col.nlx,'MarkerSize',ms);
        plot([2-dev+((b-a).*rand(1,numel(betas(:,1)))+a)],betas(:,1),'.','MarkerFaceColor', col.nacl,'MarkerEdgeColor',col.nacl,'MarkerSize',ms);
        plot([2+dev+((b-a).*rand(1,numel(betas(:,2)))+a)],betas(:,2),'.','MarkerFaceColor', col.nlx,'MarkerEdgeColor',col.nlx,'MarkerSize',ms);
        ylims = [floor(min(betas(:))-1) ceil(max(betas(:))+1)];
        %ylims = [-1.5 2.5]
        ylim([ylims]);
    end
    
    xpos = [1 2];
    y = [val_bar(1,:)-val_se(1,:); val_bar(1,:)+val_se(1,:)];
    plot([xpos;xpos], y, '-k', 'LineWidth',lw);
    plot([xpos], [val_bar(1,:) ], '.k', 'MarkerSize',ms-4);
    %errorbar([1-dev 1+dev],val_bar,val_se,'.','Color','k', 'LineWidth',lw);
    
    %set(gca,'xtick',[1 2],'xticklabel',{'Pressure','Heat'});
    ax = gca;
    ax.XAxis.FontSize = xfs;
    ax.YAxis.FontSize = afs;
    
    %legend('SAL','NLX');
     set(gca,'xtick',[1 2],'xticklabel',{'SAL','NLX'});
    title(['Coordinate: ' num2str(coordinate)]);
    ylabel('fMRI signal change (au)','FontSize',yfs);
    box off
    
    set(gcf,'color','w');
    set(gcf,'PaperPositionMode','auto');

    cd('/home/nold/Desktop/PEEP/fMRI_new/betas');
    %exportgraphics(gcf,sprintf('2betas_pharma_%d_%d_%d.png',coordinate(1),coordinate(2),coordinate(3)),'Resolution',500)
    tbl = table(betas);
    writetable(tbl,sprintf('2betas_pharma_female_%d_%d_%d.csv',coordinate(1),coordinate(2),coordinate(3)));


else
    ms = 12;
    sf = 3.9;
%     sf = 3.2;
    
    hFig = figure;
    set(hFig,'units','pixel','pos',[600 400 size1plot]);
    
    h = rm_raincloud2({betas(r.groups==1,1) betas(r.groups~=1,1); betas(r.groups==1,2) betas(r.groups~=1,2)},[col.nacl;col.nlx],'sem',ms,sf);
    
    xlim([-10 25]);
    xlabel('fMRI signal (au)','FontSize',xfs);
    set(gca,'yticklabel',{'Treatment','Baseline'},'FontSize',yfs);
   
    box off
    set(gcf,'color','w');
    set(gcf,'PaperPositionMode','auto');
    
%     print -dsvg -painters -noui /home/tinnermann/remi3/Figures/betas_raincloud_0_20_-9.svg
    print -dsvg -painters -noui /home/tinnermann/remi3/Figures/betas_raincloud_56_9_3.svg
    
end