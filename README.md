**Acute aerobic exercise does not modulate pain potentially due to differences in fitness levels and sex effects – results from a pharmacological fMRI study**
The data provided corresponds to the eLife paper (https://doi.org/10.7554/eLife.102392.1) by Nold et al. (2025). 
Description of the data and file structure

This file contains the custom behavioural and MRI analyses scripts for the study "_Acute aerobic exercise does not modulate pain potentially due to differences in fitness levels and sex effects – results from a pharmacological fMRI study_" published in _eLife_ (doi: https://doi.org/10.7554/eLife.102392.2):

- preprocsessing wrapper script (BIDS)
- first level analysis
- second level analyses
- physiology analyses/extraction
- behavioural analyses files include:
      - 01_PEEP_LOAD_IN_DATA.R # Loads in the raw data (needs to be adapted to laod in the raw data txt file only, as this also cleans the data)
      - analyses_behavioural_effects-exercise_heat.R # runs analyses (behaviour)
      - analyses_cortical_effects-exercise_heat.R # runs analyses (brain)

For a detailed, open-access preprocessing pipeline (BIDS standard) please refer to https://github.com/ChristianBuechel/spm_bids_pipeline


**The following datasets/data has been uploaded:**

 - Behavioural raw (cleaned) data

 - Con Images (uncorrected at p < 0.001) corresponding to the figures of the main article


**Files and variables**

*File: elife_behavioural_complete_data_clean.txt*

Description: 

subject: subject id

session: session id (day 2 or 3)

age: age

gender: (biological sex

modality: heat (2) or pressure (1)

pain_intensity: absolut epain intensity (in kPA or °C)

pain_rating: pain rating on VAS scale (0 - 100)

VAS: relaitve intensity of stimulus 

exercise_intensity: Exercise intensity (1 = high, 0 = low)

exercise_rating: BORG rating of exercise intensity 

exercise_block: number of exercise block (1-4)

cum_watt: cumulative watt after each block

sum_watt: sum of watts across blocks

pharm_cond: treatment condition (0 = Saline, 1 = naloxone)

nr_pain_rating: trial number (1-18)

treatment_order: order of treatment administration across days (1 = day 2 NLX, 0 = Day 2 SAL)

ftp: FTP value

pwc: weight corrected FTP value

expectation_exercise: expectation ratings on exercise 


*File: elife_figure_5_con_img_uncorrected.nii*

Description: con image (uncorrected p < 0.001) for figure 5 (interaction stimulus intensity and treatment for heat)


*File: elife_figure_2_con_img_unocrrected.nii*

Description: con image (uncorrected p < 0.001) for figure 2 (heat 70>50>30 in Saline condition for heat)


*File: elife_s_svc_mask.nii*

Description: Small volume correction mask used


*File: elife_figure_7_con_img_uncorrected.nii*

Description: con image (uncorrected p < 0.001) for figure 7 ( Exericse high intensity > exercise Low intensity in SALINE condition for heat with covariate FTP)


*File: elife_figure_8_con_img_uncorrected.nii*

Description: con image (uncorrected p < 0.001) for figure 5 (2-sample t-test between males and females for contrast interaction of exercise and drug with covariate FTP for heat)


**Code/software used**

MATLAB

SPM12 (for fMRI data)

R (packages required: svglite, reporttools, ggridges, ggdist, rstan, ggplot2, dplyr, MASS, pwr, tidyverse, rstatix, lme4, lmerTest, ggpubr, ez, reshape2, data.table, patchwork, viridis, lm.beta, car, emmeans, ggExtra, sjPlot, sjmisc, sjlabelled, merTools, ggeffects, glmmTMB, stargazer, MuMIn, plotrix, extrafont, multcomp, simr, effsize, robustlmm, plyr, lattice, readr, rmarkdown, Rmisc, tidyr, gghalves)

