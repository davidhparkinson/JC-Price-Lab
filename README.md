# JC-Price-Lab
## Instructions for Separation Code:
### If doing Peptide separation:\
1. Obtain the vials file (FRAC reg from the HPLC export), the signal file (signal from HPLC export at 214 nm, baseline 600 nm), and the baseline file (a blank run at the same wavelengths, smoothed manually to reduce any noise). Put these in their respective folders (within the HPLC folder)\
2. Run the code, indicating the # of groups and mass of peptide injected

### If doing Intact Protein Separation:\
1. Obtain the vials file (FRAC reg from the HPLC export), the signal file (signal from HPLC export at 280 nm, baseline 600 nm), and the peaks file (Integration HPLC export, at 280 nm. Look at the 3D spectrum to detect and manually unintegrated any non-proteinaceous peaks. e.g. Peaks that don't have maxima at 280 nm). Put these in their respective folders (within the HPLC folder).\
2. Run the code, indicating the # of groups and mass of protein injected. Play around with the cutoffs and width factors until you get something you like.

## Instructions for MSAnalyze Code:
1. Put all data in their respective folders in the **91Paper MS Data** folder. Coverage csv's are the proteins.csv files from PEAKS. PTM csv's are the peptide-proteins.csv files from PEAKS.
2. Run the Code and follow the input instructions.
