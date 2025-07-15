# peep_analyses_mri (MATLAB)

This file contains the custom behavioural and MRI analyses scripts for the study "_Acute aerobic exercise does not modulate pain potentially due to differences in fitness levels and sex effects – results from a pharmacological fMRI study_" published in _eLife_ (doi: https://doi.org/10.7554/eLife.102392.2):

- preprocsessing wrapper script (BIDS)
- first level analysis
- second level analyses
- physiology analyses/extraction
- behavioural analyses files include:
      - peep_compelte_data_clean.txt # cleaned raw data
      - 01_PEEP_LOAD_IN_DATA.R # Loads in the raw data (needs to be adapted to laod in the raw data txt file only, as this also cleans the data)
      - analyses_behavioural_effects-exercise_heat.R # runs analyses

For a detailed, open-access preprocessing pipeline (BIDS standard) please refer to https://github.com/ChristianBuechel/spm_bids_pipeline
