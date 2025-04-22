Matlab
# vGlut2BF_Odor_Association_Analysis

This repository contains custom MATLAB and Python scripts used in the analysis of data presented in the paper:

**"Glutamatergic projection neurons in the cholinergic basal forebrain underlie learned olfactory associational valence assignments"** (submitted to *Nature Communications*, 2025).

---

## 📁 Repository Contents
Analysis of Calcium imaging data:
- `CalciumImaging_dataanalysis/`: Scripts for processing the .csv ouput from the Inscopix software, which contains the calcium activity for each identified cell.
                                  This scripts plot the average traces, heatmap, and identify the excited, inhibited , non-responsive cells under the odor presentations during the odor screening, before and after associations.
                                  Also, it process PCA, Euclidean distance between diffe
  

- `MEDPCparser/`: Scripts for behavioral preference quantification and statistical testing
- `STDcalculator2/`: Code used to build trial-wise linear models for decoding and associative learning
- `SDTcrawl2/`: Scripts for PCA, Euclidean distance, and correlation analyses
- `dprimeRollingWindow_allmice/`: Code for generating key figures and supplementary plots
- `dprimeRollingWindow_allmice/`: Code for generating key figures and supplementary plots

---

## ⚙️ Requirements

- MATLAB R2021a or later  
- Python 3.9+  
- Required packages: `numpy`, `scipy`, `matplotlib`, `pandas`, `sklearn`

(Consider using a virtual environment or `requirements.txt` for Python dependencies.)

---

## 🚀 How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/vGlut2BF_Odor_Association_Analysis.git
