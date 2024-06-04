#=====================================================
#         LOAD IN DATA AND PREREQUISITES
#=====================================================

rm(list=ls())

library(svglite)
library(reporttools)
library(ggridges)
library(ggdist)
library(rstan)
library(ggplot2)
library(dplyr)
library(MASS)
library(pwr)
library(tidyverse)
library(rstatix)
library('lme4')
library('lmerTest')
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ez)
library(reshape2)
library(ggpubr)
library(data.table)
library(patchwork)
library(viridis)
library(lm.beta)
library(ggplot2)
library(dplyr)
library(MASS)
library(pwr)
library(tidyverse)
library(rstatix)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ez)
library(reshape2)
library(ggpubr)
library(data.table)
library(car)
library(emmeans)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(viridis)
library(lm.beta)
library(ggExtra)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(merTools)
library(ggeffects)
library(glmmTMB)
library(stargazer)
library(sjPlot)
library(sjmisc)
library(MuMIn)
library(plotrix)
library(extrafont)
library(multcomp)
library(emmeans)
library(simr)
library(effsize)
library(robustlmm)

packages <- c("plyr", "lattice", "ggplot2", "dplyr", "readr", 
              "ggplot2","rmarkdown","Rmisc", "tidyr", "gghalves")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, require, character.only = TRUE)

##############################################################
## Define Figure sizes and paths
###############################################################

save_path <- "C:/Users/user/Desktop/PEEP/fMRI/visualisation/graphs_paper_V12.0/"
supplement_path <- "C:/Users/user/Desktop/PEEP/fMRI/visualisation/graphs_supplements_V2.0/"
w_size = 10
h_size = 8
w_size_l = 12
h_size_l = 10
w_size_s = 4
h_size_s = 8
axis_text_size = 7
axis_title_size = 7
plot_title_size = 7
legend_text_size = 7
legend_title_size = 7

##############################################################
## Define path and abbreviations here 
###############################################################

path1 = 'C:/Users/user/Desktop/PEEP/fMRI/Data/LogExperiment/'
part = 'MAIN/'
subID = c('sub-01','sub-02','sub-03','sub-04','sub-06','sub-07','sub-08','sub-09','sub-11','sub-13','sub-14','sub-15','sub-17','sub-18','sub-19','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28','sub-29','sub-30',
          'sub-31','sub-32','sub-34','sub-35','sub-36','sub-38','sub-39','sub-40','sub-41','sub-43','sub-44','sub-45','sub-46','sub-47','sub-48');
csv = '.csv'
delete_first_pain_rating = FALSE # set to false if you do not want it to be deleted
ses = 'ses-01'
# Set both to FALSE when you want to look at both groups
only_control_group = FALSE # set to TRUE whether you want to only look at control group
only_experimental_group = FALSE
scale_variables = FALSE

#================================================================
# Load in Descriptives
#=================================================================
data_descriptive <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/behavioural_analysis/PEEP_descriptive_data.csv',sep = ';',header = T,dec=',')


#================================================================
# Load in Treatment List dummy coded
#=================================================================
data_treatment_dummy_coded <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/behavioural_analysis/PEEP_entblindungsliste_dummy_coded.csv',sep = ',',header = T)
data_treatment_dummy_coded<- data_treatment_dummy_coded[,-1]

# Define Treatment order: If Treatment on Day 1 allocate 2, if Treatment on Day 2 allocate 1
data_treatment_dummy_coded$treatment_order[data_treatment_dummy_coded$Day2==1] <- 1
data_treatment_dummy_coded$treatment_order[data_treatment_dummy_coded$Day2==0] <- 0
write.csv(data_treatment_dummy_coded, 'C:/Users/user/Desktop/PEEP/fMRI/Analysis/behavioural_analysis/PEEP_entblindungsliste_dummy_coded_final.csv', row.names=FALSE)
exclude_subjects <- c(5,10,12,16,20,33,37,42) # did not participate
data_treatment_dummy_coded = subset(data_treatment_dummy_coded, !(data_treatment_dummy_coded$Subject %in% exclude_subjects))


#================================================================
# Load in Exercise FTP values
#=================================================================
data_ftp <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/behavioural_analysis/ftp_values.csv',sep = ';',header = F,dec=',')
names(data_ftp)<-c("subject","FTP_value","pwc")



# ============================================================
# Data positive/negative expectation Exercise
#=============================================================
data_day1 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Data/Questionnaires/data/results-survey_day1.csv',header = TRUE,sep = ',')

# Subjects to exclude
exclude_subjects <- c("sub999","sub005","sub010","sub012","sub016","sub020","sub033","sub037","sub042")


data_exercise <- data_day1[,c('SubID',"sport.SQ001.","sport.SQ002.","sport.SQ003.","sport.SQ004.", "sport.SQ005.",
                              "sport.SQ006.","sport.SQ007.","sport.SQ008.","sport.SQ009.","sport.SQ010.","sport.SQ011.", 
                              "sport.SQ012.","sport.SQ013.","sport.SQ014." )]


# Only take items that regard pain expectation and exercise ()
data_exercise_pain <- data_exercise[,c('SubID',"sport.SQ007.","sport.SQ009.","sport.SQ014." )]
data_exercise_pain <- data_exercise_pain[complete.cases(data_exercise_pain), ]
data_exercise_pain = subset(data_exercise_pain, !(data_exercise_pain$SubID %in% exclude_subjects))
data_exercise_pain$SubID <- c(1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,38,39,40,41,43,44,45,46,47,48)
data_exercise_pain$mean_expectation <- rowMeans(data_exercise_pain[,2:4])
data_exercise_pain$expectation[data_exercise_pain$mean_expectation>4]<- 0
data_exercise_pain$expectation[data_exercise_pain$mean_expectation<4]<-  1 # put positive and neutral expectation together
data_exercise_pain$expectation[data_exercise_pain$mean_expectation==4]<- 1

#######################################################
# Load subject data exercise and pain ratings 
########################################################
all_sub <- NULL

for (i in 1:length(subID)){
  path2 <- paste(path1,part,subID[i],'/',ses,'/',sep='')
  setwd(path2)
  data <- read.csv(paste(subID[i],'_exercise_pains',csv,sep = ''),sep = ',',header = F)
  
  # make adjustments and transformations
  data$subject <- as.numeric(gsub(".*?([0-9]+).*", "\\1", subID[i]))  
  data$session <- as.factor(1)
  names(data)=c("age","gender","modality","pain_intensity","pain_rating","VAS","exercise_intensity","exercise_rating","exercise_block","cum_watt","sum_watt","pharm_cond",'subject','session')
  data$subject <- as.factor(data$subject)
  data$pain_intensity <- as.numeric(data$pain_intensity)
  data$exercise_intensity <- as.factor(data$exercise_intensity)
  data$exercise_block <- as.factor(data$exercise_block)
  data$cum_watt <- as.integer(data$cum_watt)
  data$gender <- as.factor(data$gender)
  data$nr_pain_rating <- 1:18
  data$pain_intensity <- as.factor(data$pain_intensity)
  data$treatment_order <- as.factor(data_treatment_dummy_coded$treatment_order[data_treatment_dummy_coded$Subject ==as.integer(levels(data$subject))])
  data$ftp <- data_ftp$FTP_value[data_ftp$subject ==as.integer(levels(data$subject))]
  data$pwc <- data_ftp$pwc[data_ftp$subject ==as.integer(levels(data$subject))]
  data$expectation_exercise <- as.factor(data_exercise_pain$expectation[data_exercise_pain$SubID ==as.integer(levels(data$subject))])
  data$expectation_exercise_rating <- as.factor(data_exercise_pain$mean_expectation[data_exercise_pain$SubID ==as.integer(levels(data$subject))])
  
  # assign to variable
  assign(paste0('all_', subID[i]), data)
  
  #combine
  all_sub <- rbind(all_sub,data)
}

# =============== Load session 2 Data =============================
rm(data)
all_sub_ses02 <- NULL
ses = 'ses-02'

for (i in 1:length(subID)){
  path2 <- paste(path1,part,subID[i],'/',ses,'/',sep='')
  setwd(path2)
  data <- read.csv(paste(subID[i],'_exercise_pains',csv,sep = ''),sep = ',',header = F)
  
  # make adjustments and transformations
  data$subject <- as.numeric(gsub(".*?([0-9]+).*", "\\1", subID[i]))                                   
  data$session <- as.factor(2)
  names(data)=c("age","gender","modality","pain_intensity","pain_rating","VAS","exercise_intensity","exercise_rating","exercise_block","cum_watt","sum_watt",'pharm_cond','subject','session')
  data$subject <- as.factor(data$subject)
  data$pain_intensity <- as.numeric(data$pain_intensity)
  data$exercise_intensity <- as.factor(data$exercise_intensity)
  data$exercise_block <- as.factor(data$exercise_block)
  data$cum_watt <- as.integer(data$cum_watt)
  data$gender <- as.factor(data$gender)
  data$nr_pain_rating <- 1:18
  data$pain_intensity <- as.factor(data$pain_intensity)
  data$treatment_order <- as.factor(data_treatment_dummy_coded$treatment_order[data_treatment_dummy_coded$Subject ==as.integer(levels(data$subject))])
  data$ftp <- data_ftp$FTP_value[data_ftp$subject ==as.integer(levels(data$subject))]
  data$pwc <- data_ftp$pwc[data_ftp$subject ==as.integer(levels(data$subject))]
  data$expectation_exercise <- as.factor(data_exercise_pain$expectation[data_exercise_pain$SubID ==as.integer(levels(data$subject))])
  data$expectation_exercise_rating <- as.factor(data_exercise_pain$mean_expectation[data_exercise_pain$SubID ==as.integer(levels(data$subject))])
  
  # assign to variable
  assign(paste0('all_', subID[i],'ses-02'), data)
  
  #combine
  all_sub_ses02 <- rbind(all_sub_ses02,data)
}


#===================================================================
# Combine datasets 
complete_data = rbind(all_sub,all_sub_ses02)


#===========================================================
# Clean Dataset
#==========================================================

# 1. Replace 0 with NAs (or find out rating from MRI log)
#complete_data$pain_rating %>% replace(is.na(.), 0)

# check how many NAs
sum(is.na(complete_data$pain_rating))

# check if 0s still exist
0 %in% complete_data$pain_rating
sum(complete_data$pain_rating==0,na.rm = T)

# Find out which subjects have NAs
complete_data[is.na(complete_data$pain_rating),]

# Center FTP and pwc
complete_data$ftp_c <- scale(complete_data$ftp, center = TRUE, scale = FALSE)
complete_data$pwc_c <- scale(complete_data$pwc, center = TRUE, scale = FALSE)

#============================================================
# Add new column with sporty (1) and unsporty (0)
#===============================================================
ftp_median = median(complete_data$ftp)
quartiles = summary(complete_data$ftp)

complete_data$sporty[complete_data$ftp>=ftp_median] <- 1
complete_data$sporty[complete_data$ftp<=ftp_median] <- 0

# Extreme orty /unsporty depending on quartiles
complete_data$extreme_cases_sporty[complete_data$ftp<=as.numeric(quartiles[2])] <- -0.5
complete_data$extreme_cases_sporty[complete_data$ftp>=as.numeric(quartiles[5])] <- 0.5
complete_data$extreme_cases_sporty[complete_data$ftp<=as.numeric(quartiles[5])&complete_data$ftp>=as.numeric(quartiles[2])] <- 0


# =====================================================
# Exclude sub-08 (excluded from fMRI and missing values)
#======================================================
complete_data <- complete_data[complete_data$subject!=8,]
droplevels(complete_data$subject) 



#####################################################
# ADJUST DATA FRAMES FOR LMERs AND LMMs
#####################################################

complete_data$exercise_intensity <- as.factor(complete_data$exercise_intensity)
complete_data$pharm_cond <- as.factor(complete_data$pharm_cond)
complete_data$modality <- as.factor(complete_data$modality)
complete_data$sporty <- as.factor(complete_data$sporty)
complete_data$expectation_exercise <- as.factor(complete_data$expectation_exercise)


# -------------- Heat -----------------
heat_data <- complete_data[complete_data$modality==2,]
heat_data$pain_rating_counter <- 1:9
heat_data$row <- 1:length(heat_data[,1])

# Pharm condition
heat_data_sal <- heat_data[heat_data$pharm_cond==0,]
heat_data_nlx <- heat_data[heat_data$pharm_cond==1,]

# gender
heat_data_male <- heat_data[heat_data$gender==0,]
heat_data_female <- heat_data[heat_data$gender==1,]

# -------- Pressure -------------
pressure_data <- complete_data[complete_data$modality==1,]
pressure_data$pain_rating_counter <- 1:9
pressure_data$row <- 1:length(pressure_data[,1])

# Pharma condition
pressure_data_sal <- pressure_data[pressure_data$pharm_cond==0,]
pressure_data_nlx <- pressure_data[pressure_data$pharm_cond==1,]

# Gender
pressure_data_male <- pressure_data[pressure_data$gender==0,]
pressure_data_female <- pressure_data[pressure_data$gender==1,]




##############################################################
## Define path and abbreviations here 
###############################################################

path1 = 'C:/Users/user/Desktop/PEEP/fMRI/Data/LogExperiment/'
part = 'MAIN/'
subID = c('sub-01','sub-02','sub-03','sub-04','sub-06','sub-07','sub-08','sub-09','sub-11','sub-13','sub-14','sub-15','sub-17','sub-18','sub-19','sub-21','sub-22','sub-23','sub-25','sub-26','sub-27','sub-28','sub-29','sub-30',
          'sub-31','sub-32','sub-34','sub-35','sub-36','sub-38','sub-39','sub-40','sub-41','sub-43','sub-44','sub-45','sub-46','sub-47','sub-48');
csv = '.csv'
delete_first_pain_rating = FALSE # set to false if you do not want it to be deleted
ses = 'ses-01'
# Set both to FALSE when you want to look at both groups
only_control_group = FALSE # set to TRUE whether you want to only look at control group
only_experimental_group = FALSE
scale_variables = FALSE

# --------------------
# Load in HR max data from Intervals calibration
#--------------------------------

hrmax_data <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Analysis/behavioural_analysis/HRmax.csv',sep = ',',header = F)
names(hrmax_data)<-c('subject','hrmax')

#######################################################
# LOAD DATA HR and POWER
########################################################

rm(data)
all_sub <- NULL

for (i in 1:length(subID)){
  path2 <- paste(path1,part,subID[i],'/',ses,'/',sep='')
  setwd(path2)
  data <- read.csv(paste(subID[i],'_watt_HR_exercise',csv,sep = ''),sep = ',',header = F)
  
  # make adjustments and transformations
  data$subject <- as.numeric(gsub(".*?([0-9]+).*", "\\1", subID[i]))  
  data$session <- as.factor(1)
  names(data)=c("time","watt","hr","exercise_intensity","subject","session")
  data$subject <- as.factor(data$subject)
  data$exercise_intensity <- as.factor(data$exercise_intensity)
  data$hrmax <- hrmax_data$hrmax[hrmax_data$subject ==as.integer(levels(data$subject))]
  
  # assign to variable
  assign(paste0('all_', subID[i]), data)
  
  #combine
  all_sub <- rbind(all_sub,data)
}

# =============== Load session 2 Data =============================
rm(data)
all_sub_ses02 <- NULL
ses = 'ses-02'

for (i in 1:length(subID)){
  path2 <- paste(path1,part,subID[i],'/',ses,'/',sep='')
  setwd(path2)
  data <- read.csv(paste(subID[i],'_watt_HR_exercise',csv,sep = ''),sep = ',',header = F)
  
  # make adjustments and transformations
  data$subject <- as.numeric(gsub(".*?([0-9]+).*", "\\1", subID[i]))  
  data$session <- as.factor(2)
  names(data)=c("time","watt","hr","exercise_intensity","subject","session")
  data$subject <- as.factor(data$subject)
  data$exercise_intensity <- as.factor(data$exercise_intensity)
  data$hrmax <- hrmax_data$hrmax[hrmax_data$subject ==as.integer(levels(data$subject))]
  
  # assign to variable
  assign(paste0('all_', subID[i]), data)
  
  #combine
  all_sub <- rbind(all_sub,data)
}

#===================================================================
# Combine datasets 
complete_data_hr_watt = rbind(all_sub,all_sub_ses02)

# =====================================================
# Exclude sub-08 (excluded from fMRI and missing values)
#======================================================
complete_data_hr_watt <- complete_data_hr_watt[complete_data_hr_watt$subject!=8,]
droplevels(complete_data_hr_watt$subject) 




# Defining the geom_flat_violin function ----
# Note: the below code modifies the
# existing github page by removing a parenthesis in line 50

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )



# Define a gender vector
gender_vec <- complete_data %>%
  group_by(subject,gender,pwc)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)