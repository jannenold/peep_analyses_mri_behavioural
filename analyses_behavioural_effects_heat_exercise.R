######################################################################
#
#        Exercise-induced Hypoalgesia in Heat Pain 
#                     Behavioural Analyses
#
#####################################################################

#----------------------------------------------------------------
# Aanalysis scrpit and visualisation for behavioural results
# Before running this script, the script PEEP_LOAD_IN_DATA.R 
# should be run to load in data, packages and settings.
#
# This script mainly calculates Linear Mixed Effect Models (LMER)
#
# Copyright Janne Nold 03-11-2024 (UKE, Hamburg)
#---------------------------------------------------------------


#-------------------------------------------
# Parametric Pain Ratings (Saline)
#------------------------------------------

#--------------- Statistics
lme_model_main <- lmer(pain_rating ~  VAS + treatment_order +  (1 | subject) + (1 |pain_rating_counter), data = heat_data_sal)
summary(lme_model_main)
confint(lme_model_main)

heat_data_sal$VAS <- as.factor(heat_data_sal$VAS)
lme_model_main <- lmer(pain_rating ~  VAS + treatment_order +  (1 | subject) + (1 |pain_rating_counter), data = heat_data_sal)
emmeans(lme_model_main, pairwise ~ VAS , adjust = "tukey")

# ------------- Visualisation
summary_lmer_models_main<- heat_data_sal %>%
  group_by(VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_main_sub<- heat_data_sal %>%
  group_by(subject,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(heat_data_sal, 
                          measurevar = "pain_rating", 
                          withinvars  = "VAS",
                          idvar = 'subject',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$VAS <- as.factor(summary_lmer_models_main$VAS)
summary_lmer_models_main_sub$VAS <- as.factor(summary_lmer_models_main_sub$VAS)

(param_heat_sal <-ggplot(summary_lmer_models_main,aes(VAS,pain_rating))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = VAS,y = pain_rating),shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width= 0.1),colour = 'black',fill = '#024873')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black',fill = '#024873')+
    geom_errorbar(aes(VAS,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    ylab('Pain Rating [VAS]')+xlab('Stimulus Intensity [VAS]')+ylim(0,100)+ggtitle('')+
    geom_segment(aes(x =c(0,0,0), xend = as.factor(c(30,50,70)), y = c(pain_rating), yend = c(pain_rating)),size = 1,linetype = 3, colour = "#666666")
)



#-------------------------------------------
# Interaction Treatment x Intensity 
#------------------------------------------

#----------- Statistics
lme_model_main <- lmer(pain_rating ~  VAS + pharm_cond +  treatment_order + VAS*pharm_cond + (1 | subject) + (1 | pain_rating_counter), data = heat_data)
summary(lme_model_main)
confint(lme_model_main)

#----------- Visualisation
summary_lmer_models_main<- heat_data %>%
  group_by(VAS,pharm_cond)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_main_sub<- heat_data %>%
  group_by(subject,VAS,pharm_cond)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

sum_se <- summarySEwithin(heat_data, 
                          measurevar = "pain_rating", 
                          withinvars  = c("pharm_cond","VAS"),
                          idvar = 'subject',
                          na.rm = T)


summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$VAS <- as.factor(summary_lmer_models_main$VAS)
summary_lmer_models_main_sub$VAS <- as.factor(summary_lmer_models_main_sub$VAS)

(param_heat_nlx_sal <-ggplot(summary_lmer_models_main,aes(VAS,pain_rating,fill = pharm_cond,colour = pharm_cond))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = VAS,y = pain_rating),shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.7),show.legend = F,colour = 'black')+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(VAS,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
    scale_color_manual(labels = c("SAL", "NLX"),
                       values = c("#024873", "#F29544"))+
    scale_fill_manual(labels = c("SAL", "NLX"),
                      values = c("#024873", "#F29544"))+
    theme_classic()+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('Pain Rating [VAS]')+xlab('Stimulus Intensity [VAS]')+ylim(0,100)+
    ggtitle('')+
    geom_segment(aes(x =c(0,0,0,0,0,0), xend = as.factor(c(30,30,50,50,70,70)), y = c(pain_rating), yend = c(pain_rating)),size = 0.5,linetype = 3,show.legend = F)
  
)

# Post-hoc T-tests 
heat_data$VAS <- as.factor(heat_data$VAS)
heat_data$pharm_cond <- as.factor(heat_data$pharm_cond)
lme_model_main <- lmer(pain_rating ~ VAS*pharm_cond + (1 | subject) + (1 | pain_rating_counter), data = heat_data)
emm<-emmeans(lme_model_main, pairwise ~  pharm_cond*VAS , adjust = "tukey")
eff_size(emm,sigma =sigma(lme_model_main),edf = Inf,method = "pairwise")


#---------------------------------------------------------------------
# Interaction Treatment x Intensity (Deltas NLX - SAL)
#----------------------------------------------------------------------

#---------- Statistics 

summary_lmer_models_main<- heat_data %>%
  group_by(pharm_cond,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

  summary_lmer_models_main_sd<- heat_data %>%
  group_by(pharm_cond,VAS)%>%
  summarise_at(c('pain_rating'),sd,na.rm = T)

summary_lmer_models_main_sub<- heat_data %>%
  group_by(subject,pharm_cond,VAS,treatment_order)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)


lmer_model_df <- spread(summary_lmer_models_main, pharm_cond, pain_rating)
lmer_model_df$diff_ints <- lmer_model_df$'1' -  lmer_model_df$'0'  # Low Intensity - High Intensity Data

lmer_model_df_subs <- spread(summary_lmer_models_main_sub, pharm_cond, pain_rating)
lmer_model_df_subs$diff_ints <- lmer_model_df_subs$'1' -  lmer_model_df_subs$'0'  # Low Intensity - High Intensity Data

lme_model_main <- lmer(diff_ints ~ VAS + treatment_order+ (1 | subject) , data = lmer_model_df_subs)
summary(lme_model_main)
confint(lme_model_main)

#------------- Visualisation

summary_lmer_models_main<- heat_data %>%
  group_by(pharm_cond,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_main_sub<- heat_data %>%
  group_by(subject,pharm_cond,VAS,treatment_order)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)


lmer_model_df <- spread(summary_lmer_models_main, pharm_cond, pain_rating)
lmer_model_df$diff_ints <- lmer_model_df$'1' -  lmer_model_df$'0'  # Low Intensity - High Intensity Data

lmer_model_df_subs <- spread(summary_lmer_models_main_sub, pharm_cond, pain_rating)
lmer_model_df_subs$diff_ints <- lmer_model_df_subs$'1' -  lmer_model_df_subs$'0'  # Low Intensity - High Intensity Data

sum_se <- summarySEwithin(lmer_model_df_subs, 
                          measurevar = "diff_ints", 
                          withinvars  = c("VAS"),
                          idvar = 'subject',
                          na.rm = T)

lmer_model_df$se <- sum_se$se
lmer_model_df$VAS <- as.factor(lmer_model_df$VAS)

(deltas_param <-ggplot(lmer_model_df,aes(VAS,diff_ints))+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(VAS,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
    theme_classic()+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 Pain Rating [NLX - SAL]')+xlab('Stimulus Intensity [VAS]')+
    ggtitle('')
)



# -------------------------------------------------------
# Interaction Treatment x Intensity seperated for sex
#--------------------------------------------------------

# --------------- Statistics 

# Females
lme_model_main <- lmer(pain_rating ~  VAS*pharm_cond+ treatment_order + (1 | subject)+ (1 | nr_pain_rating) , data = heat_data_female)
summary(lme_model_main)
confint(lme_model_main)

# Post-hoc T-tests 
heat_data_female$VAS <- as.factor(heat_data_female$VAS)
heat_data_female$pharm_cond <- as.factor(heat_data_female$pharm_cond)
lme_model_main <- lmer(pain_rating ~ VAS*pharm_cond + (1 | subject) + (1 | pain_rating_counter), data = heat_data_female)
emm<-emmeans(lme_model_main, pairwise ~  pharm_cond*VAS , adjust = "tukey")
eff_size(emm,sigma =sigma(lme_model_main),edf = Inf,method = "pairwise")


# Males
lme_model_main <- lmer(pain_rating ~  VAS*pharm_cond+ treatment_order + (1 | subject)+ (1 | nr_pain_rating) , data = heat_data_male)
summary(lme_model_main)
confint(lme_model_main)


# Post-hoc T-tests 
heat_data_male$VAS <- as.factor(heat_data_male$VAS)
heat_data_male$pharm_cond <- as.factor(heat_data_male$pharm_cond)
lme_model_main <- lmer(pain_rating ~ VAS*pharm_cond + (1 | subject) + (1 | pain_rating_counter), data = heat_data_male)
emm<-emmeans(lme_model_main, pairwise ~  pharm_cond*VAS , adjust = "tukey")
eff_size(emm,sigma =sigma(lme_model_main),edf = Inf,method = "pairwise")

# ---------Visualisation

summary_lmer_models<- heat_data %>%
  group_by(pharm_cond,gender,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_2<- heat_data %>%
  group_by(subject,pharm_cond,gender,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)


sum_se <- summarySEwithin(heat_data, 
                          measurevar = "pain_rating", 
                          withinvars  = c("gender","pharm_cond","VAS"),
                          idvar = 'subject',
                          na.rm = T)

summary_lmer_models$se <- sum_se$se
summary_lmer_models$pharm_cond<- as.factor(summary_lmer_models$pharm_cond)
summary_lmer_models$VAS<- as.factor(summary_lmer_models$VAS)
summary_lmer_models_2$VAS <-as.factor(summary_lmer_models_2$VAS)

# Females
(inter_mod_treat_gender_heat_female <- ggplot(summary_lmer_models[summary_lmer_models$gender==1,],aes(VAS,pain_rating,fill = pharm_cond))+
  geom_jitter(data = summary_lmer_models_2[summary_lmer_models_2$gender==1,],aes(VAS, pain_rating,fill = pharm_cond),colour = 'black',size = 0.5,alpha = 0.4,position = position_jitterdodge(dodge.width = 0.7,jitter.width = 0.1),shape = 21,show.legend = F)+
  geom_bar(stat = 'identity',position = position_dodge(0.7), width = 0.6,colour = "black",alpha = 0.8)+
  geom_errorbar(aes(VAS,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.7))+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('Pain Rating [VAS]')+xlab('Stimulus Intensity [VAS]')+
  #scale_x_discrete(limits = c("1","0"),labels = c("Female","Male"))+ggtitle('')+
  ylim(0,100)+ggtitle('Females')
)

#Males
(inter_mod_treat_gender_heat_male <- ggplot(summary_lmer_models[summary_lmer_models$gender==0,],aes(VAS,pain_rating,fill = pharm_cond))+
  geom_jitter(data =  summary_lmer_models_2[summary_lmer_models_2$gender==0,],aes(VAS, pain_rating,fill = pharm_cond),colour = 'black',size = 0.5,alpha = 0.4,position = position_jitterdodge(dodge.width = 0.7,jitter.width = 0.1),shape = 21,show.legend = F)+
  geom_bar(stat = 'identity',position = position_dodge(0.7), width = 0.6,colour = "black",alpha = 0.8)+
  geom_errorbar(aes(VAS,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width = 0.2,size = 0.5,position = position_dodge(0.7))+
  theme_classic()+
  scale_color_manual(labels = c("SAL", "NLX"),
                     values = c("#024873", "#F29544"))+
  scale_fill_manual(labels = c("SAL", "NLX"),
                    values = c("#024873", "#F29544"))+
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  guides(fill=guide_legend(title=""))+
  ylab('Pain Rating [VAS]')+xlab('Stimulus Intensity [VAS]')+
  #scale_x_discrete(limits = c("1","0"),labels = c("Female","Male"))+ggtitle('')+
  ylim(0,100)+ggtitle('Males')
)

# ---------------------------------------------------------------------------------
# Interaction Treatment x Intensity seperated for sex (Deltas (NLX -SAL))
#------------------------------------------------------------------------------------

#------ Statistics

summary_lmer_models_main<- heat_data %>%
  group_by(pharm_cond,gender,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_main_sub<- heat_data %>%
  group_by(subject,pharm_cond,gender,treatment_order,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

lmer_model_df <- spread(summary_lmer_models_main, pharm_cond, pain_rating)
lmer_model_df$diff_ints <- lmer_model_df$'1' -  lmer_model_df$'0'  # Low Intensity - High Intensity Data

lmer_model_df_subs <- spread(summary_lmer_models_main_sub, pharm_cond, pain_rating)
lmer_model_df_subs$diff_ints <- lmer_model_df_subs$'1' -  lmer_model_df_subs$'0'  # Low Intensity - High Intensity Data

lme_model_main <- lm(diff_ints ~ VAS*gender, data = lmer_model_df_subs)
summary(lme_model_main)
confint(lme_model_main)

t.test(diff_ints~gender,lmer_model_df_subs)

#--------- Visualisation

sum_se <- summarySE(lmer_model_df_subs, 
                    measurevar = "diff_ints",groupvars = c('gender'),
                    na.rm = T)

lmer_model_df$se <- sum_se$se


(deltas_gender <-ggplot(lmer_model_df,aes(gender,diff_ints,group = gender,fill =gender))+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(gender,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  = 0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    scale_color_manual(labels = c("Female", "Male"),
                       values = c("#F2F2F2", "#1A1A1A"))+
    scale_fill_manual(labels = c("Female", "Male"),
                      values = c("#F2F2F2", "#1A1A1A"))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 Pain Rating [NLX - SAL]')+xlab('')+
    ggtitle('')+
    ylim(-3,13)+
    scale_x_discrete(limits = c("1","0"),labels = c("Female","Male"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)
)



#------ Statistics

summary_lmer_models_main<- heat_data %>%
  group_by(pharm_cond,gender,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_main_sub<- heat_data %>%
  group_by(subject,pharm_cond,gender,VAS,treatment_order)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

lmer_model_df <- spread(summary_lmer_models_main, pharm_cond, pain_rating)
lmer_model_df$diff_ints <- lmer_model_df$'1' -  lmer_model_df$'0'  # Low Intensity - High Intensity Data

lmer_model_df_subs <- spread(summary_lmer_models_main_sub, pharm_cond, pain_rating)
lmer_model_df_subs$diff_ints <- lmer_model_df_subs$'1' -  lmer_model_df_subs$'0'  # Low Intensity - High Intensity Data

lme_model_main <- lmer(diff_ints ~ VAS*gender + treatment_order + (1 | subject), data = lmer_model_df_subs)
summary(lme_model_main)
confint(lme_model_main)

#--------- Visualisation

sum_se <- summarySE(lmer_model_df_subs, 
                    measurevar = "diff_ints",groupvars = c('VAS','gender'),
                    na.rm = T)

lmer_model_df$se <- sum_se$se
lmer_model_df$VAS <- as.factor(lmer_model_df$VAS)
lmer_model_df_subs$VAS <- as.factor(lmer_model_df_subs$VAS)

(deltas_gender_intensity <-ggplot(lmer_model_df,aes(VAS,diff_ints,group = gender,fill =gender))+
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(VAS,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  = 0.2,size = 0.5,position = position_dodge(0.7))+
    theme_classic()+
    scale_color_manual(labels = c("Female", "Male"),
                       values = c("#F2F2F2", "#1A1A1A"))+
    scale_fill_manual(labels = c("Female", "Male"),
                      values = c("#F2F2F2", "#1A1A1A"))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 Pain Rating [NLX - SAL]')+xlab('Stimulus Intensity [VAS]')+
    ggtitle('')+
    ylim(-3,13)+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)
)

ggsave(paste(supplement_path,'figure_deltas_param_gender.svg'), width = 6, height = 6, units = "cm")
ggsave(paste(supplement_path,'figure_deltas_param_gender.png'), width = 6, height = 6, units = "cm")

# -------------------------------------------------
# Exercise Intensity on pain ratings
#------------------------------------------------------

# On Pain Rating:
lme_model_main <- lmer(pain_rating ~   exercise_intensity + treatment_order+ (1|subject)+ (1|pain_rating_counter), data = heat_data_sal)
summary(lme_model_main)
confint(lme_model_main)

# -------------------------------------------------
# Exercise Intensity and treatment pain ratings
#------------------------------------------------------

# On Pain Rating:
lme_model_main <- lmer(pain_rating ~   exercise_intensity*pharm_cond+ treatment_order+ (1|subject)+ (1|pain_rating_counter), data = heat_data)
summary(lme_model_main)
confint(lme_model_main)

# --------------- Visualisation
summary_lmer_models_main_sub<- heat_data %>%
  group_by(subject,exercise_intensity,pharm_cond)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

summary_lmer_models_main<- heat_data %>%
  group_by(exercise_intensity,pharm_cond)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)


sum_se <- summarySEwithin(heat_data, 
                          measurevar = "pain_rating", 
                          withinvars  = c("pharm_cond","exercise_intensity"),
                          idvar = 'subject',
                          na.rm = T)

summary_lmer_models_main$se <- sum_se$se
summary_lmer_models_main$pharm_cond<- as.factor(summary_lmer_models_main$pharm_cond)
summary_lmer_models_main$VAS<- as.factor(summary_lmer_models_main$VAS)
summary_lmer_models_main_sub$pharm_cond<- as.factor(summary_lmer_models_main_sub$pharm_cond)
summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
summary_lmer_models_main_sub$VAS<- as.factor(summary_lmer_models_main_sub$VAS)

(exercise_treat_heat <-ggplot(summary_lmer_models_main,aes(exercise_intensity,pain_rating, fill = pharm_cond, colour = pharm_cond))+
    geom_jitter(data = summary_lmer_models_main_sub, aes(x = exercise_intensity,y = pain_rating),shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
    geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
    geom_errorbar(aes(exercise_intensity,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
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
    ylab('Pain Rating [VAS]')+xlab('Exercise Intensity')+ylim(0,100)+ggtitle('')+
    scale_x_discrete(limits = c("0","1"),labels = c("Low","High"))
    
   
)





# -------------------------------------------------
# Exercise Intensity and fitness level on pain ratings
#------------------------------------------------------

# On Pain Rating:
lme_model_main <- lmer(pain_rating ~   exercise_intensity*pwc+ treatment_order+ (1|subject)+ (1|pain_rating_counter), data = heat_data_sal)
summary(lme_model_main)
confint(lme_model_main)


# -------------------------------------------------
# Fitness Level on EIH 
#------------------------------------------------------

#-----------------SALINE 

# On difference pain ratings:  
summary_lmer_models_main<- heat_data_sal %>%
  group_by(subject,exercise_intensity,pwc,treatment_order,gender)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

lmer_model_df <- spread(summary_lmer_models_main, exercise_intensity, pain_rating)
lmer_model_df$diff_hi_low_rating <- lmer_model_df$'0' -  lmer_model_df$'1'  # Low Intensity - High Intensity Data


# fit regular model
lme_model_main <- lm(diff_hi_low_rating ~   pwc + treatment_order, data = lmer_model_df)
summary(lme_model_main)
confint(lme_model_main)

# Test Correlation
cor.test(lmer_model_df$diff_hi_low_rating,lmer_model_df$pwc)


#---------------- Visualisation

(saline_fitness <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = '#024873',fill = '#024873'))+
    geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
    theme_classic()+
    scale_color_manual(labels = c("SAL"),
                       values = c("#024873"))+
    scale_fill_manual(labels = c("SAL"),
                      values = c("#024873"))+
    geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
    theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
          axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
          plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),
          strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
    ggtitle('Saline')
)


# -------------------------------------------------
# Fitness Level and Sex on EIH 
#------------------------------------------------------

#------------------- SALINE

lme_model_main <- lm(diff_hi_low_rating ~ treatment_order+ pwc*gender, data = lmer_model_df)
summary(lme_model_main)
confint(lme_model_main)

cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==0],lmer_model_df$pwc[lmer_model_df$gender==0])
cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==1],lmer_model_df$pwc[lmer_model_df$gender==1])

(saline_fitness_gender <-  
    ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,fill = gender,colour = gender))+
    geom_point(size = 0.5,shape = 21,alpha = 0.9,show.legend = F)+
    geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend =T, fullrange = T)+
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
    ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Saline')
)

# -------------------------------------------------
# Fitness Level and Sex and Treatment on EIH 
#------------------------------------------------------

# ------- Statistics 

summary_lmer_models_main<- heat_data %>%
  group_by(subject,exercise_intensity,pwc,gender,treatment_order,pharm_cond,VAS)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)
lmer_model_df <- spread(summary_lmer_models_main, exercise_intensity, pain_rating)
lmer_model_df$diff_hi_low_rating <- lmer_model_df$'0' -  lmer_model_df$'1'  # Low Intensity - High Intensity Data


lme_model_main <- lmer(diff_hi_low_rating ~   VAS +treatment_order + pharm_cond*gender*pwc+(1|subject), data = lmer_model_df)
summary(lme_model_main)
confint(lme_model_main)



# fit robust model:
lme_model_main_robust <- rlmer(diff_hi_low_rating ~   VAS+treatment_order + pharm_cond*gender*pwc+(1|subject), data = lmer_model_df)
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
model_males <- lmer_model_df[lmer_model_df$gender == 0,]
model_males_saline <- model_males[model_males$pharm_cond==0,]
cor.test(model_males_saline$diff_hi_low_rating,model_males_saline$pwc)

# Females
model_females <- lmer_model_df[lmer_model_df$gender == 1,]
model_females_saline <- model_females[model_females$pharm_cond==0,]
cor.test(model_females_saline$diff_hi_low_rating,model_females_saline$pwc)

# NLX:

# Males: 
model_males <- lmer_model_df[lmer_model_df$gender == 0,]
model_males_nlx <- model_males[model_males$pharm_cond==1,]
cor.test(model_males_nlx$diff_hi_low_rating,model_males_nlx$pwc)

# Females
model_females <- lmer_model_df[lmer_model_df$gender == 1,]
model_females_nlx <- model_females[model_females$pharm_cond==1,]
cor.test(model_females_nlx$diff_hi_low_rating,model_females_nlx$pwc)

#-------------------
# Note: Labels reversed since femlaes are 1 and males are 0
# (in fMRI data Males are -1 and Females are 1)
lmer_model_df_sal<-lmer_model_df[lmer_model_df$pharm_cond==0,]

(eih_gender_saline <- 
    ggplot(lmer_model_df_sal,aes(x = pwc,y=diff_hi_low_rating,fill = gender,colour = gender))+
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
    ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Saline')
)

lmer_model_df_nlx<-lmer_model_df[lmer_model_df$pharm_cond==1,]


(eih_gender_naloxone <- 
    ggplot(lmer_model_df_nlx,aes(x = pwc,y=diff_hi_low_rating,fill = gender,colour = gender))+
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
    ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
    xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
    ggtitle('Naloxone')
)


#-------------------------------------------
# Exercise Intensity: Power (Low vs. High)
#-------------------------------------------

#---------------- Statistics 
summary_data_watt <- complete_data_hr_watt %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('watt'),mean,na.rm = T)

t.test(summary_data_watt$watt[summary_data_watt$exercise_intensity == 1],summary_data_watt$watt[summary_data_watt$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(summary_data_watt$watt[summary_data_watt$exercise_intensity == 1],summary_data_watt$watt[summary_data_watt$exercise_intensity == 0],paired=TRUE)


#------------------------Visualisation
plot_data<- summary_data_watt
plot_data_add <- spread(plot_data, exercise_intensity, watt)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
#mean_diff <- c(mean_diff,mean_diff)
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$watt[plot_data$exercise_intensity == 0]
dataRight   <-  plot_data$watt[plot_data$exercise_intensity== 1]

# Change y axis range
y_lim_min   <-  0
y_lim_max   <-  250

# Change colors for plot

leftColor         <-  "#005C53"
rightColor        <-  "#3C008E"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

power_cycling <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("Power [Watt]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(11, y_lim_max))+
  ggtitle("")

power_cycling

#-------------------------------------------
# Exercise Intensity: Heart Rate (Low vs. High)
#-------------------------------------------

# -------------- Statistics 

summary_data_hr <- complete_data_hr_watt %>%
  group_by(subject,exercise_intensity,hrmax,session)%>%
  summarise_at(c('hr'),mean,na.rm = T)

# HR did not work in subject:
# remove sub-024
summary_data_hr <- summary_data_hr[summary_data_hr$subject!=24, ]

# remove sub-036
summary_data_hr<- summary_data_hr[summary_data_hr$subject!=36, ]

# remove sub-030
summary_data_hr<- summary_data_hr[summary_data_hr$subject!=30, ]

# remove sub-032
summary_data_hr<- summary_data_hr[summary_data_hr$subject!=32, ]

# remove sub-034
summary_data_hr<- summary_data_hr[summary_data_hr$subject!=34, ]

summary_data_hr_plot <- summary_data_hr %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('hr'),mean,na.rm = T)

# paired t test
t.test(summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 1],summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 1],summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 0],paired=TRUE)


#---------------------------------------------------------

plot_data<- summary_data_hr_plot
plot_data_add <- spread(plot_data, exercise_intensity, hr)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
#mean_diff <- c(mean_diff,mean_diff)
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$hr[plot_data$exercise_intensity == 0]
dataRight   <-  plot_data$hr[plot_data$exercise_intensity== 1]

# Change y axis range
y_lim_min   <-  0
y_lim_max   <-  180

# Change colors for plot
leftColor         <-  "#005C53"
rightColor        <-  "#3C008E"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

hr_cycling <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("Heart Rate [bpm]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(11, y_lim_max))+
  ggtitle("")

hr_cycling


mean(dataLeft)
mean(dataRight)

#-------------------------------------------
# Exercise Intensity: RPE (Borg Scale) (Low vs. High)
#-------------------------------------------

# --------- Statistics 

exercise_data<- complete_data[,c('subject','pharm_cond','exercise_intensity','exercise_rating','sporty')]

# set ratings below 7 to 7
exercise_data$exercise_rating[exercise_data$exercise_rating < 7] <- 7


exercise_data_plot<- exercise_data %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('exercise_rating'),mean,na.rm = T)

# By Intensity
low_data <- exercise_data[exercise_data$exercise_intensity==0,]
high_data <- exercise_data[exercise_data$exercise_intensity==1,]

# paired t test
t.test(exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 1],exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 1],exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 0],paired=TRUE)


# ---------- RAINCLOUD PLOT ------------

plot_data<- exercise_data_plot
plot_data_add <- spread(plot_data, exercise_intensity, exercise_rating)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$exercise_rating[plot_data$exercise_intensity == 0]
dataRight   <-  plot_data$exercise_rating[plot_data$exercise_intensity== 1]

# Change y axis range
y_lim_min   <-  6
y_lim_max   <-  20

# Change colors for plot
leftColor         <-  "#005C53"
rightColor        <-  "#3C008E"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

borg_rating <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("Rating of Perceived Exertion\n[BORG 6-20]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(6, y_lim_max))+
  ggtitle("")

borg_rating



################################################################################################
######       COMBINE FIGURES TO PUBLICATION FIGURES 
################################################################################################

#---------------- NEW FIGURE: Pain Ratings EIH with Treatment
ggarrange(exercise_treat_heat)

ggsave(paste(save_path,'figure_extra.svg'), width = 7, height = 7, units = "cm")
ggsave(paste(save_path,'figure_extra.png'), width = 7, height = 7, units = "cm")

# -------------- Figure 2: Sanity Checks Pain  
ggarrange(ggarrange(param_heat_sal, ncol = 1, nrow = 1,labels = c("A"),common.legend = T,legend = 'right',font.label = list(size = 11),align = 'hv',widths = c(0.55,1.25)),
          ncol = 1,nrow = 1,heights = c(1,0.75))

ggsave(paste(save_path,'figure_2.svg'), width = 6, height = 8, units = "cm")
ggsave(paste(save_path,'figure_2.png'), width = 6, height = 8, units = "cm")


# -------------- Figure 3: Sanity Checks Exercise  
ggarrange(power_cycling,hr_cycling, borg_rating,ncol = 3, nrow = 1,labels = c("A","B","C"),font.label = list(size = 11),widths = c(0.75,0.75,0.75))

ggsave(paste(save_path,'figure_3.svg'), width = 15, height = 6, units = "cm")
ggsave(paste(save_path,'figure_3.png'), width = 15, height = 6, units = "cm")

# -------------- Figure 4: Treatent effect Naloxone
ggarrange(ggarrange(param_heat_nlx_sal,ggarrange(inter_mod_treat_gender_heat_female,inter_mod_treat_gender_heat_male+rremove('ylab'),common.legend = T,legend = 'right', nrow = 1,ncol = 2,labels = c("C",'D'),font.label = list(size = 11)),ncol = 1, nrow = 2,labels = c("A"),font.label = list(size = 11)),
          ggarrange(deltas_param,deltas_gender_intensity+rremove('legend'),labels = c("B","E"),font.label = list(size = 11),ncol = 1,nrow = 2,align='hv'),
          ncol = 2,nrow= 1,widths = c(1.25,0.75))

ggsave(paste(save_path,'figure_4.svg'), width = 14, height = 13, units = "cm")
ggsave(paste(save_path,'figure_4.png'), width = 14, height = 13, units = "cm")


# -------------- Figure 6: Exercise induced Hypoalgesia rACC

ggarrange(saline_fitness,common.legend = T,legend = 'none',labels = c("A"),font.label = list(size = 11),align = 'hv')

ggsave(paste(save_path,'figure_6a.svg'), width = 4.5, height = 6, units = "cm")
ggsave(paste(save_path,'figure_6a.png'), width = 4.5, height = 6, units = "cm")

# -------------- Figure 7: Sex x Pharma x Fitness

ggarrange(eih_gender_saline,eih_gender_naloxone,common.legend = T,legend = 'bottom',labels = c("A","B"),font.label = list(size = 11),align = 'hv',ncol = 1)

ggsave(paste(save_path,'figure_7.svg'), width = 4.5, height = 10, units = "cm")
ggsave(paste(save_path,'figure_7.png'), width = 4.5, height = 10, units = "cm")

