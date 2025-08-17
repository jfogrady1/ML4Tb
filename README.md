# Accurate and robust classification of _Mycobacterium bovis_-infected cattle using peripheral blood RNA-seq data

---

## Summary
The zoonotic bacterium, _Mycobacterium bovis_, causes bovine tuberculosis (bTB) and is closely related to _Mycobacterium tuberculosis_, the primary cause of human tuberculosis (hTB). Bovine TB remains recalcitrant to eradication in endemic countries where current diagnostics fail to identify all infected animals. While blood-based RNA biomarkers identified through machine learning have shown accurate discrimination of hTB-positive and hTB-negative individuals, similar approaches have not been explored for bTB. Here, we use RNA-seq and machine learning to investigate the utility of peripheral blood mRNA as a host-response biomarker for bTB using data from Ireland, the UK and the US. We identify a 13-gene signature and a 273-gene elastic net classifier that differentiate bTB-positive from bTB-negative cattle, achieving area under the curve (AUC) values of 0.994/0.902 for the former and 0.967/0.942 for the latter in training and testing, respectively. Both classifiers demonstrate high sensitivity (≥0.912) and specificity (0.853) in the testing set. Additionally, we show that they robustly distinguish bTB+ animals from those infected with other bacterial or viral pathogens (AUC ≥0.821).  These RNA-based classifiers accurately diagnose bTB and differentiate bTB from other diseases, representing a promising tool for augmenting current diagnostics to advance bTB eradication efforts in endemic regions.

---

## Overview of Study

<img src="https://github.com/jfogrady1/ML4Tb/blob/45f5d25d76a54c9f8641ed4fbe3807f976dd1d67/Figure_01.png" alt ="Overview">

---

## Data availability

| Study | Abbreviation | Infectious agent/<br>Disease | Accession Number(s) | 
| ----- | ----- | ----- | ----- | 
| O'Grady et al. (2025) | OGR25-BTB | *M. bovis* (bTB) | GSE255724 | 
| McLoughlin et al. (2014) | MCL14-BTB | *M. bovis* (bTB) | PRJNA257841 | 
| McLoughlin et al. (2021) | MCL21-BTB | *M. bovis* (bTB) | PRJEB27764<br>PRJEB44470 | 
| Wiarda et al. (2020) | WIA20-BTB | *M. bovis* (bTB) | PRJNA600004 | 
| Alonso-Hearn et al. (2019) | ALO19-MAP | *M. avium ssp. paratuberculosis* (Johne's disease) | PRJNA565369 | 
| Johnston et al. (2021) | JOH21-BRD | BRSV (BRD) | GSE152959 | 
| O'Donoghue et al. (2023) | ODO23-BRD | BoHV-1 (BRD) | GSE199108 | 

- Imputed and filtered WGS data from O'Grady et al. (2025) is available at Zenodo ((https://zenodo.org/records/13752056). 
- Raw genotype information derived from the RNA-seq variant call for bTB studies described here is available at Zenodo (https://doi.org/10.5281/zenodo.16887104). 