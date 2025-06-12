# 생물정보학 및 실습1 - Bioinformatics and Lab Practice 1

This repository contains assignments and practice notebooks for the SNU course **생물정보학 및 실습1**.  
The main notebook in this repository demonstrates data analysis, including feature counting and basic preprocessing.

---

## 📁 Contents

- `CoLab_TermProj_2025_1.ipynb` — The first week's Jupyter notebook for the lab work
- `CoLab_TermProj_2025_2.ipynb` — The second week's Jupyter notebook for the lab work
- `CoLab_TermProj_2025_3.ipynb` — The third week's Jupyter notebook for the lab work. The resulting figures for third week can be found in `fig/week3`
- `MOA_prompt.ipynb` & `MOA.R`  — Final term project notebook
- `analysis_results/` — Final figures and analysis outputs

---

## 💻 Plans for Own Analysis
- Analyze the results of CLIP performed with a protein other than LIN28 using the same method. When applying Guided Mission 3, I may need to adjust the methods depending on the experimental approach or protein type.

### Progress
- As an initial step, I applied the entropy and base mutation quantification pipeline from Guided Mission 3 to a CLIP-seq dataset of IGF2BP1 focusing on the ACTB gene. This served as a prototype to test my ability to reuse core methods on a different RBP. I will upload final analysis code and ppt next week.

# Final Term Project My Own Analysis: Translational Profiling of LIN28A-bound Transcripts

## Goal
To investigate how **LIN28A binding influences ribosome occupancy** around start codons, reusing methods from Guided Mission 3.

---

## Datasets used
- Ribo-seq datasets: 
	-`fivepcounts-filtered-RPF-siLuc.txt` (control) 
	-`fivepcounts-filtered-RPF-siLin28a.txt` (LIN28A knockdown)
- Exon annotations: Start and stop codon regions (Gencode)
- CLIP targets (from LIN28A): Used to stratify transcripts

---

## Methods Summary
1. Filtered 5′ read ends aligned to start codons ±150 nt  
2. Binned and normalized read counts by transcript  
3. Compared average ribosome density across conditions (±SEM)  
4. Per-position t-tests with FDR correction  
5. Created heatmaps for transcript-level occupancy  
6. Quantified pause score changes and visualized volcano-style plot  

All core pipelines were built by modifying and expanding the code from the guided missions.

---

## Key Figure
- 📈 **Mean Ribosome Density Plot** (Figure S5A style):  
  Shows loss of start codon ribosome pausing after LIN28A knockdown

- 🔥 **Transcript-Level Heatmap**:  
  Heatmap showing reduced pausing across LIN28A-regulated transcripts

- 🌋 **Pause Score Volcano Plot**:  
  Highlights transcripts with significant fold-change in pausing behavior

All plots are saved in `analysis_results/figures/`.

---
### Final Deliverables

- `MOA_prompt.ipynb` & `MOA.R` — Full analysis pipeline
- `analysis_results/figures/` — Final figures (PNG)
- `README.md` — This file

---

📅 Project completed: **June 2025**  
👩🏻‍💻 Author: Yerim
