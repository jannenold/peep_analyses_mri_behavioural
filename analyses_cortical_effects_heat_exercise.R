######################################################################
#
#        Exercise-induced Hypoalgesia in Heat Pain 
#                     Cortical Analyses
#
#####################################################################

#----------------------------------------------------------------
# Aanalysis scrpt and visualisation for cortical results
# Before running this script, the script PEEP_LOAD_IN_DATA.R 
# should be run to load in data and settings
#
#
# Copyright Janne Nold 03-11-2024 (UKE, Hamburg)
#---------------------------------------------------------------


#----------------------------------------------------------
#   Parametric Effect Saline Heat: betas
# ---------------------------------------------------------

# -------- Insula MNI: 36 6 14

betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/6betas_intensity_heat_int_sal_36_6_14.csv',sep = ',',header = T)
names(betas) <- c('h30_sal','h30_nlx','h50_sal','h50_nlx','h70_sal','h70_nlx')
betas <- betas[,c('h30_sal','h50_sal',"h70_sal")]
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]
betas$gender <- gender_vec$gender
betas_long <- gather(betas,group,beta,h30_sal:h70_sal)
betas_long$VAS[betas_long$group == 'h30_sal'] <- 30
betas_long$VAS[betas_long$group == 'h50_sal'] <- 50
betas_long$VAS[betas_long$group == 'h70_sal'] <- 70


betas_long$VAS <- as.factor(betas_long$VAS)
summary_lmer_models_main<- betas_long %>%
  group_by(VAS)%>%
  summarise_at(c('beta'),mean,na.rm = T)

summary_lmer_models_main_sub<- betas_long %>%
  group_by(SubID,VAS)%>%
  summarise_at(c('beta'),mean,na.rm = T)

# Calculating SE for within design (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = "VAS",
                          idvar = 'SubID',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$VAS <- as.factor(summary_lmer_models_main$VAS)
summary_lmer_models_main_sub$VAS <- as.factor(summary_lmer_models_main_sub$VAS)

(param_heat_sal_insula <-ggplot(summary_lmer_models_main,aes(VAS,beta))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = VAS,y = beta),shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width= 0),colour = 'black',fill = '#024873')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black',fill = '#024873')+
    geom_errorbar(aes(VAS,ymin = beta-se,ymax=beta+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('fMRI signal change [au]')+xlab('Stimulus Intensity [VAS]')+ggtitle('')+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)
)

#----- dpIns MNI: 39 -15 18

betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/6betas_intensity_heat_int_sal_39_-15_18.csv',sep = ',',header = T)
names(betas) <- c('h30_sal','h50_sal','h70_sal')
betas <- betas[,c('h30_sal','h50_sal',"h70_sal")]
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]
betas$gender <- gender_vec$gender
betas_long <- gather(betas,group,beta,h30_sal:h70_sal)
betas_long$VAS[betas_long$group == 'h30_sal'] <- 30
betas_long$VAS[betas_long$group == 'h50_sal'] <- 50
betas_long$VAS[betas_long$group == 'h70_sal'] <- 70


betas_long$VAS <- as.factor(betas_long$VAS)
summary_lmer_models_main<- betas_long %>%
  group_by(VAS)%>%
  summarise_at(c('beta'),mean,na.rm = T)

summary_lmer_models_main_sub<- betas_long %>%
  group_by(SubID,VAS)%>%
  summarise_at(c('beta'),mean,na.rm = T)

# Calculating SE for within design (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = "VAS",
                          idvar = 'SubID',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$VAS <- as.factor(summary_lmer_models_main$VAS)
summary_lmer_models_main_sub$VAS <- as.factor(summary_lmer_models_main_sub$VAS)

(param_heat_sal_dpIns <-ggplot(summary_lmer_models_main,aes(VAS,beta))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = VAS,y = beta),shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width= 0),colour = 'black',fill = '#024873')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black',fill = '#024873')+
    geom_errorbar(aes(VAS,ymin = beta-se,ymax=beta+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('fMRI signal change [au]')+xlab('Stimulus Intensity [VAS]')+ggtitle('')+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)
)



#  ----------- middle cingulate cortex MNI: 6, 10, 39

betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/6betas_intensity_heat_int_sal_6_10_39.csv',sep = ',',header = T)
names(betas) <- c('h30_sal','h50_sal','h70_sal')
betas <- betas[,c('h30_sal','h50_sal',"h70_sal")]
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]
betas$gender <- gender_vec$gender
betas_long <- gather(betas,group,beta,h30_sal:h70_sal)
betas_long$VAS[betas_long$group == 'h30_sal'] <- 30
betas_long$VAS[betas_long$group == 'h50_sal'] <- 50
betas_long$VAS[betas_long$group == 'h70_sal'] <- 70
betas_long$VAS <- as.factor(betas_long$VAS)


summary_lmer_models_main<- betas_long %>%
  group_by(VAS)%>%
  summarise_at(c('beta'),mean,na.rm = T)

summary_lmer_models_main_sub<- betas_long %>%
  group_by(SubID,VAS)%>%
  summarise_at(c('beta'),mean,na.rm = T)

# Calculating SE for within design (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = "VAS",
                          idvar = 'SubID',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$VAS <- as.factor(summary_lmer_models_main$VAS)
summary_lmer_models_main_sub$VAS <- as.factor(summary_lmer_models_main_sub$VAS)

(param_heat_sal_mcc <-ggplot(summary_lmer_models_main,aes(VAS,beta))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = VAS,y = beta),shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width= 0),colour = 'black',fill = '#024873')+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',fill = '#024873')+
    geom_errorbar(aes(VAS,ymin = beta-se,ymax=beta+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('fMRI signal change [au]')+xlab('Stimulus Intensity [VAS]')+ggtitle('')+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)
)

# -----------------------------------------
# Treatment Effect: HEAT NLX > HEAT SAL
# -----------------------------------------

# ---- Insula 46 8 6
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/6betas_intensity_heat_sal_nlx_46_8_6.csv',sep = ',',header = T)
names(betas) <- c('betas_sal_h30','betas_nlx_h30','betas_sal_h50','betas_nlx_h50','betas_sal_h70','betas_nlx_h70')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, betas_sal_h30:betas_nlx_h70)

betas_long$intensity[betas_long$group == 'betas_sal_h30'] <- 30
betas_long$intensity[betas_long$group == 'betas_nlx_h30'] <- 30
betas_long$intensity[betas_long$group == 'betas_sal_h50'] <- 50
betas_long$intensity[betas_long$group == 'betas_nlx_h50'] <- 50
betas_long$intensity[betas_long$group == 'betas_sal_h70'] <- 70
betas_long$intensity[betas_long$group == 'betas_nlx_h70'] <- 70

betas_long$pharm[betas_long$group == 'betas_nlx_h30'] <- 1
betas_long$pharm[betas_long$group == 'betas_sal_h30'] <- 0
betas_long$pharm[betas_long$group == 'betas_nlx_h50'] <- 1
betas_long$pharm[betas_long$group == 'betas_sal_h50'] <- 0
betas_long$pharm[betas_long$group == 'betas_nlx_h70'] <- 1
betas_long$pharm[betas_long$group == 'betas_sal_h70'] <- 0

betas_long$intensity <- as.factor(betas_long$intensity)
betas_long$pharm <- as.factor(betas_long$pharm)

summary_betas<- betas_long %>%
  group_by(intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)

sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)


summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$intensity <- as.factor(summary_betas$intensity)


(heat_nlx_sal_insula <- ggplot(summary_betas,aes(intensity,beta,fill = pharm,colour = pharm))+
    geom_jitter(data= betas_long,aes(x = intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.4),
                show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
    geom_errorbar(aes(intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.4))+
    geom_point(data= summary_betas,aes(x = intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
    theme_classic()+
    scale_color_manual(labels = c("SAL", "NLX"),
                       values = c("#024873", "#F29544"))+
    scale_fill_manual(labels = c("SAL", "NLX"),
                      values = c("#024873", "#F29544"))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('fMRI signal change [au]')+xlab('')+
    scale_x_discrete(limits = c("30","50","70"),labels = c("30","50","70"))+
    xlab('Stimulus Intensity [VAS]')+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    ylim(-4.5,6.25)+
    ggtitle(''))


#----- Deltas between NLX - SAL

summary_betas<- betas_long %>%
  group_by(intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)

summary_betas_sub<- betas_long %>%
  group_by(SubID,intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)

lmer_model_df <- spread(summary_betas, pharm, beta)
lmer_model_df$diff_ints <- lmer_model_df$'1' -  lmer_model_df$'0'  # Low Intensity - High Intensity Data

lmer_model_df_subs <- spread(summary_betas_sub, pharm, beta)
lmer_model_df_subs$diff_ints <- lmer_model_df_subs$'1' -  lmer_model_df_subs$'0'  # Low Intensity - High Intensity Data

sum_se <- summarySEwithin(lmer_model_df_subs, 
                          measurevar = "diff_ints", 
                          withinvars  = c("intensity"),
                          idvar = 'SubID',
                          na.rm = T)

lmer_model_df$se <- sum_se$se


(deltas_param_antIns <-ggplot(lmer_model_df,aes(intensity,diff_ints))+
    # geom_point(data = lmer_model_df_subs, aes(x = VAS,y = diff_ints),shape = 21,alpha = 0.4,size = 0.5,show.legend = F,colour = 'black')+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(intensity,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
    theme_classic()+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 fMRI signal [NLX - SAL]')+xlab('Stimulus Intensity [VAS]')+
    #ylim(0,8)+
    ggtitle('')
)



# -------FIR model
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/FIR_heat_nlx_sal_insula.csv',sep = ',',header = T)
names(betas) <- c('time','h70_sal','h70_nlx','h70_sal_se','h70_nlx_se')

betas_long <- gather(betas,group,beta,h70_sal:h70_nlx)

betas_long$pharm[betas_long$group == 'h70_nlx'] <- 1
betas_long$pharm[betas_long$group == 'h70_nlx_se'] <- 1
betas_long$pharm[betas_long$group == 'h70_sal'] <- 0
betas_long$pharm[betas_long$group == 'h70_sal_se'] <- 0

betas_long_2 <- gather(betas_long,group2,se,h70_sal_se:h70_nlx_se)
betas_long_2$pharm <- as.factor(betas_long_2$pharm)

heat_insula_fir <-ggplot(betas_long_2,aes(time,beta,fill = pharm,colour = pharm))+
  annotate('rect', xmin=5, xmax=18, ymin=min(betas_long_2$beta), ymax=max(betas_long_2$beta)+0.1, alpha=.2, fill='#666666')+
  geom_ribbon(aes(ymin = beta-se,ymax = beta+se),alpha = 0.3,show.legend = F,size = 0.001)+
  geom_line(size = 1,show.legend = F)+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873","#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c( "#024873","#F29544"))+
  
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm")) +
  guides(fill=guide_legend(title="")) +
  #ylim(-0.13,0.4)+
  xlim(-0.5,20)+
  ylab('fMRI signal change [au]')+xlab('Time [sec]')+
  #geom_vline(xintercept = 0,colour = '#666666',size = 1)+
  geom_vline(xintercept = 0,colour = '#666666',size = 0.5)+
  geom_vline(xintercept = 18,colour = '#666666',size = 0.5)+
  annotate("text", x = c(2,15),y = c(max(betas_long_2$beta)+0.075,max(betas_long_2$beta)+0.075),label = c("Stimulus\nStart","Stimulus\nEnd"),size = 2)+ggtitle("")


heat_insula_fir


# ------------------------------------
# Interaction Intensity x Treatment Heat
# -------------------------------------

# PAG -2 -24 -8
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/6betas_intensity_treat_-2_-24_-8.csv',sep = ',',header = T)
names(betas) <- c('betas_sal_h30','betas_nlx_h30','betas_sal_h50','betas_nlx_h50','betas_sal_h70','betas_nlx_h70')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, betas_sal_h30:betas_nlx_h70)

betas_long$intensity[betas_long$group == 'betas_sal_h30'] <- 30
betas_long$intensity[betas_long$group == 'betas_nlx_h30'] <- 30
betas_long$intensity[betas_long$group == 'betas_sal_h50'] <- 50
betas_long$intensity[betas_long$group == 'betas_nlx_h50'] <- 50
betas_long$intensity[betas_long$group == 'betas_sal_h70'] <- 70
betas_long$intensity[betas_long$group == 'betas_nlx_h70'] <- 70


betas_long$pharm[betas_long$group == 'betas_nlx_h30'] <- 1
betas_long$pharm[betas_long$group == 'betas_sal_h30'] <- 0
betas_long$pharm[betas_long$group == 'betas_nlx_h50'] <- 1
betas_long$pharm[betas_long$group == 'betas_sal_h50'] <- 0
betas_long$pharm[betas_long$group == 'betas_nlx_h70'] <- 1
betas_long$pharm[betas_long$group == 'betas_sal_h70'] <- 0

betas_long$intensity <- as.factor(betas_long$intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)


summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$intensity <- as.factor(summary_betas$intensity)



inter_treat_int_heat_pag <- ggplot(summary_betas,aes(intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.4),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.4))+
  geom_point(data= summary_betas,aes(x = intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  #geom_bar(stat = 'identity',position = position_dodge(0.7), width = 0.6,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("30","50","70"),labels = c("30","50","70"))+
  xlab('Stimulus Intensity [VAS]')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  ylim(-4.5,6.25)+
  ggtitle('')

inter_treat_int_heat_pag

#----- Deltas between NLX - SAL

summary_betas<- betas_long %>%
  group_by(intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)

summary_betas_sub<- betas_long %>%
  group_by(SubID,intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)

lmer_model_df <- spread(summary_betas, pharm, beta)
lmer_model_df$diff_ints <- lmer_model_df$'1' -  lmer_model_df$'0'  # Low Intensity - High Intensity Data

lmer_model_df_subs <- spread(summary_betas_sub, pharm, beta)
lmer_model_df_subs$diff_ints <- lmer_model_df_subs$'1' -  lmer_model_df_subs$'0'  # Low Intensity - High Intensity Data

sum_se <- summarySEwithin(lmer_model_df_subs, 
                          measurevar = "diff_ints", 
                          withinvars  = c("intensity"),
                          idvar = 'SubID',
                          na.rm = T)

lmer_model_df$se <- sum_se$se


(deltas_param_PAG <-ggplot(lmer_model_df,aes(intensity,diff_ints))+
    # geom_point(data = lmer_model_df_subs, aes(x = VAS,y = diff_ints),shape = 21,alpha = 0.4,size = 0.5,show.legend = F,colour = 'black')+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(intensity,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
    theme_classic()+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 fMRI signal [NLX - SAL]')+xlab('Stimulus Intensity [VAS]')+
    #ylim(0,8)+
    ggtitle('')
)


# ----- FIR model
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/FIR_inter_treat_int_heat_pag_-2_-24_-8.csv',sep = ',',header = T)
names(betas) <- c('time','h70_sal','h70_nlx','h70_sal_se','h70_nlx_se')

betas_long <- gather(betas,group,beta,h70_sal:h70_nlx)

betas_long$pharm[betas_long$group == 'h70_nlx'] <- 1
betas_long$pharm[betas_long$group == 'h70_nlx_se'] <- 1
betas_long$pharm[betas_long$group == 'h70_sal'] <- 0
betas_long$pharm[betas_long$group == 'h70_sal_se'] <- 0

betas_long_2 <- gather(betas_long,group2,se,h70_sal_se:h70_nlx_se)
betas_long_2$pharm <- as.factor(betas_long_2$pharm)

inter_treat_int_heat_pag_fir <-ggplot(betas_long_2,aes(time,beta,fill = pharm,colour = pharm))+
  annotate('rect', xmin=5, xmax=18, ymin=min(betas_long_2$beta), ymax=max(betas_long_2$beta)+0.1, alpha=.2, fill='#666666')+
  geom_ribbon(aes(ymin = beta-se,ymax = beta+se),alpha = 0.3,show.legend = F,size = 0.001)+
  geom_line(size = 1,show.legend = F)+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873","#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c( "#024873","#F29544"))+
  
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm")) +
  guides(fill=guide_legend(title="")) +
  #ylim(-0.13,0.4)+
  xlim(-0.5,20)+
  ylab('fMRI signal change [au]')+xlab('Time [sec]')+
  #geom_vline(xintercept = 0,colour = '#666666',size = 1)+
  geom_vline(xintercept = 0,colour = '#666666',size = 0.5)+
  geom_vline(xintercept = 18,colour = '#666666',size = 0.5)+
  annotate("text", x = c(2,15),y = c(max(betas_long_2$beta)+0.075,max(betas_long_2$beta)+0.075),label = c("Stimulus\nStart","Stimulus\nEnd"),size = 2)+ggtitle("")


inter_treat_int_heat_pag_fir


# ------------------------------------
# Exercise Intensity on Pain (SAL): Exercise high > low
# -------------------------------------

# -8 20 -30
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_-8_20_-30.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')


# 8 -8 33
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_8_-8_33.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')

# ------------------------------------
# Exercise Intensity on Pain (SAL): Exercis low > high (sal)
# -------------------------------------

# 30 -15 15
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_30_-15_15.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')


# 30 8 -14
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_30_8_-14.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')


#--------------------------------------
# ROI Analysis
#-------------------------------------

# ---- PAG (Brainstem NAV)
dataexLo30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exLo_brainstemnav_PAG.csv',sep = ',',header = T)
dataexLo50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exLo_brainstemnav_PAG.csv',sep = ',',header = T)
dataexLo70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exLo_brainstemnav_PAG.csv',sep = ',',header = T)

dataexHi30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exHi_brainstemnav_PAG.csv',sep = ',',header = T)
dataexHi50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exHi_brainstemnav_PAG.csv',sep = ',',header = T)
dataexHi70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exHi_brainstemnav_PAG.csv',sep = ',',header = T)

data_comb <- cbind(dataexLo30,dataexLo50[,7],dataexLo70[,7],dataexHi30[,7],dataexHi50[,7],dataexHi70[,7])


# SAL only
data_comb_sal <- data_comb[data_comb$pharm == -1,]
data_comb_sal$subid <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
data_comb_sal$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

data_comb_sal <- data_comb_sal[,7:14]
names(data_comb_sal) <- c("lo30","lo50","lo70","hi30","hi50","hi70","sub","treatment_order")
datacomb_sal_long <- gather(data_comb_sal,group,betas,lo30:hi70)

datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo30'] <- 30
datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo50'] <- 50
datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo70'] <- 70

datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi30'] <- 30
datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi50'] <- 50
datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi70'] <- 70

datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo30'] <- 0
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo50'] <- 0
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo70'] <- 0

datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi30'] <- 1
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi50'] <- 1
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi70'] <- 1

lme_model_main <- lmer(betas ~   exercise_intensity + treatment_order+ (1|sub), data = datacomb_sal_long)
summary(lme_model_main)
confint(lme_model_main)

# Both Treatment conditions
dataexLo30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exLo_brainstemnav_PAG.csv',sep = ',',header = T)
dataexLo50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exLo_brainstemnav_PAG.csv',sep = ',',header = T)
dataexLo70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exLo_brainstemnav_PAG.csv',sep = ',',header = T)

dataexHi30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exHi_brainstemnav_PAG.csv',sep = ',',header = T)
dataexHi50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exHi_brainstemnav_PAG.csv',sep = ',',header = T)
dataexHi70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exHi_brainstemnav_PAG.csv',sep = ',',header = T)

data_comb <- cbind(dataexLo30,dataexLo50[,7],dataexLo70[,7],dataexHi30[,7],dataexHi50[,7],dataexHi70[,7])
data_comb$subid <- rep(c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48),2)
data_comb$treatment_order <- rep(data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8],2)

data_comb <- data_comb[,c(1,7:14)]
names(data_comb) <- c("pharm","lo30","lo50","lo70","hi30","hi50","hi70","sub","treatment_order")
data_comb_long <- gather(data_comb,group,betas,lo30:hi70)

data_comb_long$VAS[data_comb_long$group == 'lo30'] <- 30
data_comb_long$VAS[data_comb_long$group == 'lo50'] <- 50
data_comb_long$VAS[data_comb_long$group == 'lo70'] <- 70

data_comb_long$VAS[data_comb_long$group == 'hi30'] <- 30
data_comb_long$VAS[data_comb_long$group == 'hi50'] <- 50
data_comb_long$VAS[data_comb_long$group == 'hi70'] <- 70

data_comb_long$exercise_intensity[data_comb_long$group == 'lo30'] <- 0
data_comb_long$exercise_intensity[data_comb_long$group == 'lo50'] <- 0
data_comb_long$exercise_intensity[data_comb_long$group == 'lo70'] <- 0

data_comb_long$exercise_intensity[data_comb_long$group == 'hi30'] <- 1
data_comb_long$exercise_intensity[data_comb_long$group == 'hi50'] <- 1
data_comb_long$exercise_intensity[data_comb_long$group == 'hi70'] <- 1


lme_model_main <- lmer(betas ~   exercise_intensity*pharm + treatment_order+ (1|sub), data = data_comb_long)
summary(lme_model_main)
confint(lme_model_main)


# --------------- Visualisation
summary_lmer_models_main_sub<- data_comb_long %>%
  group_by(sub,exercise_intensity,pharm)%>%
  summarise_at(c('betas'),mean,na.rm = T)

summary_lmer_models_main<- data_comb_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('betas'),mean,na.rm = T)


sum_se <- summarySEwithin(data_comb_long, 
                          measurevar = "betas", 
                          withinvars  = c("pharm","exercise_intensity"),
                          idvar = 'sub',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$pharm<- as.factor(summary_lmer_models_main$pharm)
summary_lmer_models_main_sub$pharm<- as.factor(summary_lmer_models_main_sub$pharm)
summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)

(exercise_treat_heat_pag <-ggplot(summary_lmer_models_main,aes(exercise_intensity,betas, fill = pharm, colour = pharm))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = exercise_intensity,y = betas),shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
    geom_errorbar(aes(exercise_intensity,ymin = betas-se,ymax=betas+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    scale_color_manual(labels = c("SAL", "NLX"),
                       values = c("#024873", "#F29544"))+
    scale_fill_manual(labels = c("SAL", "NLX"),
                      values = c("#024873", "#F29544"))+
    guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('fMRI signal change [au]')+xlab('Exercise Intensity')+ggtitle('')+
    scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
    geom_hline(yintercept = 0,width = 0.5)
  
)

# ---- RVM (Tinnermann)
dataexLo30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exLo_RVM.csv',sep = ',',header = T)
dataexLo50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exLo_RVM.csv',sep = ',',header = T)
dataexLo70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exLo_RVM.csv',sep = ',',header = T)

dataexHi30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exHi_RVM.csv',sep = ',',header = T)
dataexHi50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exHi_RVM.csv',sep = ',',header = T)
dataexHi70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exHi_RVM.csv',sep = ',',header = T)

data_comb <- cbind(dataexLo30,dataexLo50[,7],dataexLo70[,7],dataexHi30[,7],dataexHi50[,7],dataexHi70[,7])


# SAL only
data_comb_sal <- data_comb[data_comb$pharm == -1,]
data_comb_sal$subid <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
data_comb_sal$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

data_comb_sal <- data_comb_sal[,7:14]
names(data_comb_sal) <- c("lo30","lo50","lo70","hi30","hi50","hi70","sub","treatment_order")
datacomb_sal_long <- gather(data_comb_sal,group,betas,lo30:hi70)

datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo30'] <- 30
datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo50'] <- 50
datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo70'] <- 70

datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi30'] <- 30
datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi50'] <- 50
datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi70'] <- 70

datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo30'] <- 0
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo50'] <- 0
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo70'] <- 0

datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi30'] <- 1
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi50'] <- 1
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi70'] <- 1

lme_model_main <- lmer(betas ~   exercise_intensity + treatment_order+ (1|sub), data = datacomb_sal_long)
summary(lme_model_main)
confint(lme_model_main)

# Both Treatment conditions
dataexLo30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exLo_RVM.csv',sep = ',',header = T)
dataexLo50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exLo_RVM.csv',sep = ',',header = T)
dataexLo70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exLo_RVM.csv',sep = ',',header = T)

dataexHi30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exHi_RVM.csv',sep = ',',header = T)
dataexHi50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exHi_RVM.csv',sep = ',',header = T)
dataexHi70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exHi_RVM.csv',sep = ',',header = T)

data_comb <- cbind(dataexLo30,dataexLo50[,7],dataexLo70[,7],dataexHi30[,7],dataexHi50[,7],dataexHi70[,7])
data_comb$subid <- rep(c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48),2)
data_comb$treatment_order <- rep(data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8],2)

data_comb <- data_comb[,c(1,7:14)]
names(data_comb) <- c("pharm","lo30","lo50","lo70","hi30","hi50","hi70","sub","treatment_order")
data_comb_long <- gather(data_comb,group,betas,lo30:hi70)

data_comb_long$VAS[data_comb_long$group == 'lo30'] <- 30
data_comb_long$VAS[data_comb_long$group == 'lo50'] <- 50
data_comb_long$VAS[data_comb_long$group == 'lo70'] <- 70

data_comb_long$VAS[data_comb_long$group == 'hi30'] <- 30
data_comb_long$VAS[data_comb_long$group == 'hi50'] <- 50
data_comb_long$VAS[data_comb_long$group == 'hi70'] <- 70

data_comb_long$exercise_intensity[data_comb_long$group == 'lo30'] <- 0
data_comb_long$exercise_intensity[data_comb_long$group == 'lo50'] <- 0
data_comb_long$exercise_intensity[data_comb_long$group == 'lo70'] <- 0

data_comb_long$exercise_intensity[data_comb_long$group == 'hi30'] <- 1
data_comb_long$exercise_intensity[data_comb_long$group == 'hi50'] <- 1
data_comb_long$exercise_intensity[data_comb_long$group == 'hi70'] <- 1


lme_model_main <- lmer(betas ~   exercise_intensity*pharm + treatment_order+ (1|sub), data = data_comb_long)
summary(lme_model_main)
confint(lme_model_main)


# --------------- Visualisation
summary_lmer_models_main_sub<- data_comb_long %>%
  group_by(sub,exercise_intensity,pharm)%>%
  summarise_at(c('betas'),mean,na.rm = T)

summary_lmer_models_main<- data_comb_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('betas'),mean,na.rm = T)


sum_se <- summarySEwithin(data_comb_long, 
                          measurevar = "betas", 
                          withinvars  = c("pharm","exercise_intensity"),
                          idvar = 'sub',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$pharm<- as.factor(summary_lmer_models_main$pharm)
summary_lmer_models_main_sub$pharm<- as.factor(summary_lmer_models_main_sub$pharm)
summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)

(exercise_treat_heat_rvm <-ggplot(summary_lmer_models_main,aes(exercise_intensity,betas, fill = pharm, colour = pharm))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = exercise_intensity,y = betas),shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
    geom_errorbar(aes(exercise_intensity,ymin = betas-se,ymax=betas+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    scale_color_manual(labels = c("SAL", "NLX"),
                       values = c("#024873", "#F29544"))+
    scale_fill_manual(labels = c("SAL", "NLX"),
                      values = c("#024873", "#F29544"))+
    guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('fMRI signal change [au]')+xlab('Exercise Intensity')+ggtitle('')+
    scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
    geom_hline(yintercept = 0,width = 0.5)
  
)

# ---- frontal midline
dataexLo30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exLo_frontMidline_vega.csv',sep = ',',header = T)
dataexLo50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exLo_frontMidline_vega.csv',sep = ',',header = T)
dataexLo70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exLo_frontMidline_vega.csv',sep = ',',header = T)

dataexHi30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exHi_frontMidline_vega.csv',sep = ',',header = T)
dataexHi50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exHi_frontMidline_vega.csv',sep = ',',header = T)
dataexHi70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exHi_frontMidline_vega.csv',sep = ',',header = T)

data_comb <- cbind(dataexLo30,dataexLo50[,7],dataexLo70[,7],dataexHi30[,7],dataexHi50[,7],dataexHi70[,7])


# SAL only
data_comb_sal <- data_comb[data_comb$pharm == -1,]
data_comb_sal$subid <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
data_comb_sal$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

data_comb_sal <- data_comb_sal[,7:14]
names(data_comb_sal) <- c("lo30","lo50","lo70","hi30","hi50","hi70","sub","treatment_order")
datacomb_sal_long <- gather(data_comb_sal,group,betas,lo30:hi70)

datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo30'] <- 30
datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo50'] <- 50
datacomb_sal_long$VAS[datacomb_sal_long$group == 'lo70'] <- 70

datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi30'] <- 30
datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi50'] <- 50
datacomb_sal_long$VAS[datacomb_sal_long$group == 'hi70'] <- 70

datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo30'] <- 0
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo50'] <- 0
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'lo70'] <- 0

datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi30'] <- 1
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi50'] <- 1
datacomb_sal_long$exercise_intensity[datacomb_sal_long$group == 'hi70'] <- 1

lme_model_main <- lmer(betas ~   exercise_intensity + treatment_order+ (1|sub), data = datacomb_sal_long)
summary(lme_model_main)
confint(lme_model_main)

# Both Treatment conditions
dataexLo30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exLo_frontMidline_vega.csv',sep = ',',header = T)
dataexLo50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exLo_frontMidline_vega.csv',sep = ',',header = T)
dataexLo70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exLo_frontMidline_vega.csv',sep = ',',header = T)

dataexHi30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat30_exHi_frontMidline_vega.csv',sep = ',',header = T)
dataexHi50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat50_exHi_frontMidline_vega.csv',sep = ',',header = T)
dataexHi70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_rois_analysis/roi_betas_heat70_exHi_frontMidline_vega.csv',sep = ',',header = T)


data_comb <- cbind(dataexLo30,dataexLo50[,7],dataexLo70[,7],dataexHi30[,7],dataexHi50[,7],dataexHi70[,7])
data_comb$subid <- rep(c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48),2)
data_comb$treatment_order <- rep(data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8],2)

data_comb <- data_comb[,c(1,7:14)]
names(data_comb) <- c("pharm","lo30","lo50","lo70","hi30","hi50","hi70","sub","treatment_order")
data_comb_long <- gather(data_comb,group,betas,lo30:hi70)

data_comb_long$VAS[data_comb_long$group == 'lo30'] <- 30
data_comb_long$VAS[data_comb_long$group == 'lo50'] <- 50
data_comb_long$VAS[data_comb_long$group == 'lo70'] <- 70

data_comb_long$VAS[data_comb_long$group == 'hi30'] <- 30
data_comb_long$VAS[data_comb_long$group == 'hi50'] <- 50
data_comb_long$VAS[data_comb_long$group == 'hi70'] <- 70

data_comb_long$exercise_intensity[data_comb_long$group == 'lo30'] <- 0
data_comb_long$exercise_intensity[data_comb_long$group == 'lo50'] <- 0
data_comb_long$exercise_intensity[data_comb_long$group == 'lo70'] <- 0

data_comb_long$exercise_intensity[data_comb_long$group == 'hi30'] <- 1
data_comb_long$exercise_intensity[data_comb_long$group == 'hi50'] <- 1
data_comb_long$exercise_intensity[data_comb_long$group == 'hi70'] <- 1


lme_model_main <- lmer(betas ~   exercise_intensity*pharm + treatment_order+ (1|sub), data = data_comb_long)
summary(lme_model_main)
confint(lme_model_main)


# --------------- Visualisation
summary_lmer_models_main_sub<- data_comb_long %>%
  group_by(sub,exercise_intensity,pharm)%>%
  summarise_at(c('betas'),mean,na.rm = T)

summary_lmer_models_main<- data_comb_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('betas'),mean,na.rm = T)


sum_se <- summarySEwithin(data_comb_long, 
                          measurevar = "betas", 
                          withinvars  = c("pharm","exercise_intensity"),
                          idvar = 'sub',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$pharm<- as.factor(summary_lmer_models_main$pharm)
summary_lmer_models_main_sub$pharm<- as.factor(summary_lmer_models_main_sub$pharm)
summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)

(exercise_treat_heat_frontmid <-ggplot(summary_lmer_models_main,aes(exercise_intensity,betas, fill = pharm, colour = pharm))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = exercise_intensity,y = betas),shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
    geom_errorbar(aes(exercise_intensity,ymin = betas-se,ymax=betas+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    scale_color_manual(labels = c("SAL", "NLX"),
                       values = c("#024873", "#F29544"))+
    scale_fill_manual(labels = c("SAL", "NLX"),
                      values = c("#024873", "#F29544"))+
    guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('fMRI signal change [au]')+xlab('Exercise Intensity')+ggtitle('')+
    scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
    geom_hline(yintercept = 0,width = 0.5)
  
)

# ------------------------------------
# Exercise Intensity and treatment on Pain: Interaction
# -------------------------------------

# -4 -33 -27
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_-4_-33_-27.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')


# 20 -36 -42
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_20_-36_-42.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')

 #------------------------------------
  # Exercise Intensity and treatment on Pain: Interaction
  # -------------------------------------

# 9 -8 32
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_9_-8_32.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')


# 16 39 -16
betas <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/4betas_exercise_treatment_16_39_-16.csv',sep = ',',header = T)
names(betas) <- c('heat_lo_sal','heat_lo_nlx','heat_hi_sal','heat_hi_nlx')
betas$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
betas$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_long <- gather(betas,group,beta, heat_lo_sal:heat_hi_nlx)

betas_long$exercise_intensity[betas_long$group == 'heat_lo_sal'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$exercise_intensity[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$exercise_intensity[betas_long$group == 'heat_hi_sal'] <- 1


betas_long$pharm[betas_long$group == 'heat_lo_sal'] <- 1
betas_long$pharm[betas_long$group == 'heat_lo_nlx'] <- 0
betas_long$pharm[betas_long$group == 'heat_hi_nlx'] <- 1
betas_long$pharm[betas_long$group == 'heat_hi_sal'] <- 0

betas_long$exercise_intensity <- as.factor(betas_long$exercise_intensity)
betas_long$pharm <- as.factor(betas_long$pharm)


summary_betas<- betas_long %>%
  group_by(exercise_intensity,pharm)%>%
  summarise_at(c('beta'),mean,na.rm = T)


sum_se <- summarySEwithin(betas_long, 
                          measurevar = "beta", 
                          withinvars  = c("exercise_intensity","pharm"),
                          idvar = 'SubID',
                          na.rm = T)

summary_betas$se <- sum_se$se
summary_betas$pharm <- as.factor(summary_betas$pharm)
summary_betas$exercise_intensity <- as.factor(summary_betas$exercise_intensity)


ggplot(summary_betas,aes(exercise_intensity,beta,fill = pharm,colour = pharm))+
  #geom_line(data= summary_betas,aes(x = intensity,y = beta,group = pharm),show.legend = F,alpha = 0.5,size = 0.75,colour = 'darkgrey',position = position_dodge(0.4))+
  geom_jitter(data= betas_long,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm,group = pharm),shape = 21,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.5),
              show.legend = F,alpha = 0.2,size = 0.3,shape = 21)+
  geom_errorbar(aes(exercise_intensity,ymin = beta-se,ymax = beta+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.5))+
  #geom_point(data= summary_betas,aes(x = exercise_intensity,y = beta,fill = pharm,colour = pharm),shape = 21,position = position_dodge(0.4),alpha = 1,size = 0.5,colour = "black")+
  geom_bar(stat = 'identity',position = position_dodge(0.5), width = 0.4,colour = "black",alpha = 0.8)+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('fMRI signal change [au]')+xlab('')+
  scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))+
  xlab('Exercise Intensity')+
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #ylim(-4.5,6.25)+
  ggtitle('')


#-----------------------------------------------------------------------
# Fitness Level on EIH (SAL): rACC MNI: 6 45  10
#-----------------------------------------------------------------------

data <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/cor_betas_exercise_heat_sal_6_45_10.csv',sep = ',',header = T)
names(data) <- c('low_ex','hi_exercise')
data$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
data$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]
data$gender <- gender_vec$gender
data$ftp  <- gender_vec$pwc
data_long <- gather(data,exercise_intensity,beta,low_ex:hi_exercise)

#----- Visualisation absolute betas
(exercise_intensity_ftp_rACC <- ggplot(data_long,aes(ftp,beta,fill = exercise_intensity,colour = exercise_intensity))+
  geom_point(size = 0.5,alpha  =0.5,show.legend = T)+
  geom_smooth(method = 'lm',alpha = 0.2,size = 1,se = T,show.legend = F)+
  scale_color_manual(labels = c("High ",'Low'),
                     values = c("#3C008E", "#005C53"))+
  scale_fill_manual(labels = c("High",'Low'),
                    values = c("#3C008E", "#005C53"))+
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title="")) +
  #ylim(-2,5.2)+
  ylab('fMRI signal change [au]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+ggtitle("")
)


# --------- Calculate deltas between betas of LI and Hi exercise

data$diff_hi_low_betas <- data$low_ex -  data$hi_exercise  # Low Intensity - High Intensity Data

(exercise_intensity_ftp_rACC_delta <- ggplot(data,aes(ftp,diff_hi_low_betas,fill = '#024873',colour = '#024873'))+
  geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
  geom_smooth(method = 'lm',alpha = 0.2,size = 1,se = T,show.legend = F)+
  theme_classic()+
  scale_color_manual(labels = c("SAL"),
                     values = c("#024873"))+
  scale_fill_manual(labels = c("SAL"),
                    values = c("#024873"))+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title="")) +
  geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
  ylab('\u0394 fMRI signal change [LI - HI exercise]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+ggtitle("Saline")

)

#-----------------------------------------------------------------------
# Fitness Level and Sex on EIH (SAL): rACC MNI: 6 45  10
#-----------------------------------------------------------------------

(exercise_intensity_ftp_rACC_delta_gender<-  ggplot(data,aes(ftp,diff_hi_low_betas,fill = gender,colour = gender))+
    geom_point(size = 0.5,shape = 21,alpha = 0.9,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F, fullrange = T)+
    theme_classic()+
    scale_color_manual(labels = c("Female", "Male"),
                       values = c("#1C02C7", "#C75302"))+
    scale_fill_manual(labels = c("Female", "Male"),
                      values = c("#1C02C7", "#C75302"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('\u0394 fMRI signal change [LI - HI exercise]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Saline')
  
)



#-----------------------------------------------------------------------
# Fitness Level and Sex and Pharma on EIH: 
#-----------------------------------------------------------------------

# ------ vmPFV MNI: 12 64 2

data <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/fmri_analysis/betas/8betas_inter_treat_exercise_gender_12_64_2.csv',sep = ',',header = T)
names(data) <- c('low_ex_sal','hi_ex_sal','low_ex_nlx','hi_ex_nlx')
data$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
data$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]
data$gender <- gender_vec$gender
data$ftp  <- gender_vec$pwc
data$diff_hi_low_betas_sal <- data$low_ex_sal -  data$hi_ex_sal  # Low Intensity - High Intensity Data
data$diff_hi_low_betas_nlx <- data$low_ex_nlx -  data$hi_ex_nlx  # Low Intensity - High Intensity Data
data_long <- gather(data,pharm_cond,diff_betas,diff_hi_low_betas_sal:diff_hi_low_betas_nlx)


#---- Visualisation 
lmer_model_df_sal<-data_long[data_long$pharm_cond=='diff_hi_low_betas_sal',]

(eih_gender_saline_vmPFC <- 
    ggplot(lmer_model_df_sal,aes(x = ftp,y=diff_betas,fill = gender,colour = gender))+
    geom_point(size = 0.5,shape = 21,alpha = 0.9,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = T, fullrange = T)+
    theme_classic()+
    scale_color_manual(labels = c("Female", "Male"),
                       values = c("#1C02C7", "#C75302"))+
    scale_fill_manual(labels = c("Female", "Male"),
                      values = c("#1C02C7", "#C75302"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    #guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('\u0394 fMRI signal change\n[LI - HI Exercise]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Saline')
)

lmer_model_df_nlx<-data_long[data_long$pharm_cond=='diff_hi_low_betas_nlx',]


(eih_gender_naloxone_vmPFC <- 
    ggplot(lmer_model_df_nlx,aes(x = ftp,y=diff_betas,fill = gender,colour = gender))+
    geom_point(size = 0.5,shape = 21,alpha = 0.9,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = T, fullrange = T)+
    theme_classic()+
    scale_color_manual(labels = c("Female", "Male"),
                       values = c("#1C02C7", "#C75302"))+
    scale_fill_manual(labels = c("Female", "Male"),
                      values = c("#1C02C7", "#C75302"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    #guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('\u0394 fMRI signal change\n[LI - HI Exercise]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Naloxone')
)


#----------------------------------------------
# ---- PAG Eippert ROI 0 -35 -17 5mm
#-------------------------------------------

dataexLo_30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_ROIs/roi_betas_heat30_exLo_PAG_0_-35_-17_5mm.csv',sep = ',',header = T)
dataexLo_50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_ROIs/roi_betas_heat50_exLo_PAG_0_-35_-17_5mm.csv',sep = ',',header = T)
dataexLo_70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_ROIs/roi_betas_heat70_exLo_PAG_0_-35_-17_5mm.csv',sep = ',',header = T)

dataexHi_30 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_ROIs/roi_betas_heat30_exHi_PAG_0_-35_-17_5mm.csv',sep = ',',header = T)
dataexHi_50 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_ROIs/roi_betas_heat50_exHi_PAG_0_-35_-17_5mm.csv',sep = ',',header = T)
dataexHi_70 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/ROI/betas_ROIs/roi_betas_heat70_exHi_PAG_0_-35_-17_5mm.csv',sep = ',',header = T)

dataexLo <- rbind(dataexLo_30,dataexLo_50,dataexLo_70)
dataexHi <- rbind(dataexHi_30,dataexHi_50,dataexHi_70)

names(dataexLo) <- c("pharm","ftp","gender","modality","int","ex_int","betas_lo")
names(dataexHi) <- c("pharm","ftp","gender","modality","int","ex_int","betas_hi")

dataexLo$gender <- as.factor(dataexLo$gender)
dataexLo$pharm <- as.factor(dataexLo$pharm)
dataexLo$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
dataexLo$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

dataexHi$gender <- as.factor(dataexHi$gender)
dataexHi$pharm <- as.factor(dataexHi$pharm)
dataexHi$SubID <- c(1,2,3,4,6,7,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
dataexHi$treatment_order <- data_treatment_dummy_coded$treatment_order[!data_treatment_dummy_coded$Subject==8]

betas_hi <- dataexHi$betas_hi
data_new<- cbind(dataexLo,betas_hi)
data_new$diff_hi_low_betas <- data_new$betas_lo -  data_new$betas_hi  # Low Intensity - High Intensity Data

#- -------- Statistics

# LMER
lme_model_main <- lmer(diff_hi_low_betas ~   int + ftp*gender*pharm+ treatment_order+ (1|SubID), data = data_new)
summary(lme_model_main)
confint(lme_model_main)
write.csv(round(summary(lme_model_main)$coefficients,2),'C:/Users/user/Desktop/xxx.csv')

# Robust LMER
lme_model_main_robust <- rlmer(diff_hi_low_betas ~   int +treatment_order + gender*pharm*ftp+(1|SubID), data = data_new)
summary(lme_model_main_robust)

# get coefficients from non-robust model to extract Satterthwaite approximated DFs
coefs <- data.frame(coef(summary(lme_model_main)))

# get coefficients from robust model to extract t-values
coefs.robust <- coef(summary(lme_model_main_robust))

# calculate p-values based on robust t-values and non-robust approx. DFs
p.values <- 2*pt(abs(coefs.robust[,3]), coefs$df, lower=FALSE)
p.values


# Correlation

# SALINE:
# Males: 
model_males <- data_new[data_new$gender == -1,]
model_males_saline <- model_males[model_males$pharm==-1,]
cor.test(model_males_saline$diff_hi_low_betas,model_males_saline$ftp)

# Females
model_females <- data_new[data_new$gender == 1,]
model_females_saline <- model_females[model_females$pharm==-1,]
cor.test(model_females_saline$diff_hi_low_betas,model_females_saline$ftp)


# NLX:
# Males: 
# Males: 
model_males <- data_new[data_new$gender == -1,]
model_males_nlx <- model_males[model_males$pharm==1,]
cor.test(model_males_nlx$diff_hi_low_betas,model_males_nlx$ftp)

# Females
model_females <- data_new[data_new$gender == 1,]
model_females_nlx <- model_females[model_females$pharm==1,]
cor.test(model_females_nlx$diff_hi_low_betas,model_females_nlx$ftp)

#------------ Visualisation

lmer_model_df_sal<-data_new[data_new$pharm==-1,]

(eih_gender_saline_roi_pag <- 
    ggplot(lmer_model_df_sal,aes(x = ftp,y=diff_hi_low_betas,fill = gender,colour = gender))+
    geom_point(size = 0.5,shape = 21,alpha = 0.9,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = T, fullrange = T)+
    theme_classic()+
    scale_color_manual(labels = c("Male", "Female"),
                       values = c("#C75302", "#1C02C7"))+
    scale_fill_manual(labels = c("Male", "Female"),
                      values = c("#C75302", "#1C02C7"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    #guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('\u0394 fMRI signal change\n[LI - HI Exercise]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Saline')
)

lmer_model_df_nlx<-data_new[data_new$pharm==1,]


(eih_gender_naloxone_roi_pag <- 
    ggplot(lmer_model_df_nlx,aes(x = ftp,y=diff_hi_low_betas,fill = gender,colour = gender))+
    geom_point(size = 0.5,shape = 21,alpha = 0.9,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = T, fullrange = T)+
    theme_classic()+
    scale_color_manual(labels = c("Male", "Female"),
                       values = c("#C75302", "#1C02C7"))+
    scale_fill_manual(labels = c("Male", "Female"),
                      values = c("#C75302", "#1C02C7"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    #guides(fill=guide_legend(title=""))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = legend_title_size),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('\u0394 fMRI signal change\n[LI - HI Exercise]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Naloxone')
)

ggarrange(eih_gender_saline,eih_gender_naloxone,common.legend = T,legend = 'right')



################################################################################################
######       COMBINE FIGURES TO PUBLICATION FIGURES 
################################################################################################


# ----------- Figure 2
ggarrange(param_heat_sal_insula,param_heat_sal_dpIns,param_heat_sal_mcc, ncol = 1, nrow = 3,labels = c("E","F","G"),common.legend = T,legend = 'right',font.label = list(size = 11),align = 'hv',widths = c(0.75,1.25))

ggsave(paste(save_path,'figure_2a.svg'), width = 4, height = 13, units = "cm")
ggsave(paste(save_path,'figure_2a.png'), width = 4, height = 13, units = "cm")


#----------- Figure 5

ggarrange(ggarrange(inter_treat_int_heat_pag,deltas_param_PAG,inter_treat_int_heat_pag_fir,common.legend = T,legend = 'none',widths = c(0.85,0.8,1.5),labels = c("B","C","D"),font.label = list(size = 11),align = 'hv',ncol = 3),
         # ggarrange(heat_nlx_sal_insula,deltas_param_antIns,heat_insula_fir,common.legend = T,legend = 'none',widths = c(0.85,0.8,1.5),labels = c("f","g","h"),font.label = list(size = 11),align = 'hv',ncol = 3),
          nrow = 1,common.legend = T,legend = 'none',align = 'hv')

ggsave(paste(save_path,'figure_5.svg'), width = 14, height = 4.5, units = "cm")
ggsave(paste(save_path,'figure_5.png'), width = 14, height = 4.5, units = "cm")


#----------- Figure 5 Supplements

ggarrange(ggarrange(heat_nlx_sal_insula,deltas_param_antIns,heat_insula_fir,common.legend = T,legend = 'none',widths = c(0.85,0.8,1.5),labels = c("B","C","D"),font.label = list(size = 11),align = 'hv',ncol = 3),
             nrow = 1,common.legend = T,legend = 'none',align = 'hv')

ggsave(paste(supplement_path,'figure_5_sup.svg'), width = 14, height = 4.5, units = "cm")
ggsave(paste(supplement_path,'figure_5_sup.png'), width = 14, height = 4.5, units = "cm")


#---------- Figure ROI extra
ggarrange(ggarrange(exercise_treat_heat_rvm,exercise_treat_heat_pag,exercise_treat_heat_frontmid,common.legend = T,legend = 'none',font.label = list(size = 11),align = 'hv',ncol = 1),
          nrow = 1,common.legend = T,legend = 'none',align = 'hv')

ggsave(paste(save_path,'figure_extra_roi.svg'), width = 4, height = 13, units = "cm")
ggsave(paste(save_path,'figure_extra_roi.png'), width = 4, height = 13, units = "cm")

#--------- Figure 6c

ggarrange(exercise_intensity_ftp_rACC_delta,common.legend = T,legend = 'none',labels = c("C"),font.label = list(size = 11),align = 'hv')

ggsave(paste(save_path,'figure_6c.svg'), width = 4.25, height = 6, units = "cm")
ggsave(paste(save_path,'figure_6c.png'), width = 4.25, height = 6, units = "cm")

#--------- Figure 7

ggarrange(eih_gender_saline_vmPFC,eih_gender_naloxone_vmPFC,common.legend = T,legend = 'none',labels = c("D","E"),font.label = list(size = 11),align = 'hv',ncol = 1)

ggsave(paste(save_path,'figure_7_vmPFC.svg'), width = 4.5, height = 10, units = "cm")
ggsave(paste(save_path,'figure_7_vmPFC.png'), width = 4.5, height = 10, units = "cm")

ggarrange(eih_gender_saline_roi_pag,eih_gender_naloxone_roi_pag,common.legend = T,legend = 'none',labels = c("e","f"),font.label = list(size = 11),align = 'hv',ncol = 1)

ggsave(paste(save_path,'figure_8_pag.svg'), width = 4.5, height = 10, units = "cm")
ggsave(paste(save_path,'figure_8_pag.png'), width = 4.5, height = 10, units = "cm")
