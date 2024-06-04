function plot_4betas_exercise

subs          = [1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48];
exclude    = [8];

if ~isempty(exclude)
    subs = subs(~ismember(subs,exclude));
end

base_dir = '/projects/crunchie/nold/PEEP/fMRI/Data/MAIN/rawdata/derivatives/spm_firstlevel/HRF_complete_onset_stimulus_6mm';

barplot = 1;
%plot individual data points
ind_dots = 0;
%figure size publication or not
pub = 1;

coordinate = [39 -24 36];

% Heat SAL Ex Low > ex High
% neg Cor FTP
% 62 -30 -15
% -69 -23 -9
% 6 45 9
%42 -27 15

% Heat NLX Ex Low > ex High
% neg Cor FTP
% 2 -14 5
% -16 63 -4
% -8 32 -8

%pos cor
% 39 -24 36


image_nr = [21 22]; %Heat NLX 21 22 % Heat SAL 19 20

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
ftp_median = 122;
pwc_median = 1.755;

col.base = [27 90 129]/255;
col.infs = [243 159 86]/255;

col.nacl = [255 174 73]/255;
col.remi = [22 128 164]/255;

bw = 0.8;
dev = 0.14;
err = 1; % ci 90%: 1.645

if pub
    lw = 0.9;
    size1plot = [165 135]; %for figures and publication
    afs = 10;
    xfs = 15;
    yfs = 15;
else
    lw = 0.9;
    ms = 15;
    size1plot = [320 250];
    afs = 10;
    xfs = 11;
    yfs = 12;
end

if barplot == 1
    ms = 10;
    size1plot = [185*2 135*2];
    afs = 10;
    xfs = 15;
    yfs = 15;
    sr1 = sqrt(length(ftp_vec(pwc_vec>pwc_median)));
    sr2 = sqrt(length(pwc_vec(pwc_vec<pwc_median)));

    val_bar = [mean(betas(pwc_vec>pwc_median,1)) mean(betas(pwc_vec>pwc_median,2));...
        mean(betas(pwc_vec<pwc_median,1)) mean(betas(pwc_vec<pwc_median,2))];
    val_se = [std(betas(pwc_vec>pwc_median,1))/sr1 * err std(betas(pwc_vec>pwc_median,2))/sr1 * err;...
        std(betas(pwc_vec<pwc_median,1))/sr2 * err std(betas(pwc_vec<pwc_median,2))/sr2 * err];

    hFig = figure;
    set(hFig,'units','pixel','pos',[900 400 size1plot]);

    h = bar(val_bar,'EdgeColor','none','BarWidth',bw); hold on;


    set(h(1),'FaceColor',col.base);
    set(h(2),'FaceColor',col.infs);

    if ind_dots
        set(h,'FaceAlpha',0.6)
        a = -0.05;
        b = 0.05;
        plot([1-dev+((b-a).*rand(1,numel(betas(pwc_vec>pwc_median,1)))+a)],betas(pwc_vec>pwc_median,1),'.','MarkerFaceColor', col.base,'MarkerEdgeColor',col.base,'MarkerSize',ms);
        plot([1+dev+((b-a).*rand(1,numel(betas(pwc_vec>pwc_median,2)))+a)],betas(pwc_vec>pwc_median,2),'.','MarkerFaceColor', col.infs,'MarkerEdgeColor',col.infs,'MarkerSize',ms);
        plot([2-dev+((b-a).*rand(1,numel(betas(pwc_vec<pwc_median,1)))+a)],betas(pwc_vec<pwc_median,1),'.','MarkerFaceColor', col.base,'MarkerEdgeColor',col.base,'MarkerSize',ms);
        plot([2+dev+((b-a).*rand(1,numel(betas(pwc_vec<pwc_median,2)))+a)],betas(pwc_vec<pwc_median,2),'.','MarkerFaceColor', col.infs,'MarkerEdgeColor',col.infs,'MarkerSize',ms);
        ylims = [floor(min(betas(:))-1) ceil(max(betas(:))+1)];
        ylim([ylims]);
    end

    xpos = [1-dev 1+dev 2-dev 2+dev];
    y = [val_bar(1,:)-val_se(1,:) val_bar(2,:)-val_se(2,:); val_bar(1,:)+val_se(1,:) val_bar(2,:)+val_se(2,:)];
    plot([xpos; xpos], y, '-k', 'LineWidth',lw);
    plot([xpos; xpos], [val_bar(1,:) val_bar(2,:)], '.k', 'MarkerSize',ms-4);
    % errorbar([1-dev 1+dev; 2-dev 2+dev; 3-dev 3+dev],val_bar,val_se,'.','Color','k', 'LineWidth',lw);

    set(gca,'xtick',[1 2],'xticklabel',{'Fit','Unfit'});
    ax = gca;
    ax.XAxis.FontSize = xfs;
    ax.YAxis.FontSize = afs;

    ylabel('fMRI signal change (au)','FontSize',yfs);
    box off

    legend('Exercise Low','Exercise High');

    title(['Coordinate: ' num2str(coordinate)]);


    set(gcf,'color','w');
    set(gcf,'PaperPositionMode','auto');

    cd('/home/nold/Desktop/PEEP/fMRI/Figures/');
    exportgraphics(gcf,sprintf('4betas_exercise_%d_%d_%d.png',coordinate(1),coordinate(2),coordinate(3)),'Resolution',500)

else
    ms = 12;
    sf = 3.9;
    %     sf = 3.2;

    hFig = figure;
    set(hFig,'units','pixel','pos',[600 400 size1plot]);

    h = rm_raincloud2({betas(r.groups==1,1) betas(r.groups~=1,1); betas(r.groups==1,2) betas(r.groups~=1,2)},[col.nacl;col.remi],'sem',ms,sf);

    xlim([-10 25]);
    xlabel('fMRI signal (au)','FontSize',xfs);
    set(gca,'yticklabel',{'Treatment','Baseline'},'FontSize',yfs);

    box off
    set(gcf,'color','w');
    set(gcf,'PaperPositionMode','auto');

    %     print -dsvg -painters -noui /home/tinnermann/remi3/Figures/betas_raincloud_0_20_-9.svg
    print -dsvg -painters -noui /home/tinnermann/remi3/Figures/betas_raincloud_56_9_3.svg

end