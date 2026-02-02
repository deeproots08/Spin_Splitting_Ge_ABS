# Spin_Splitting_Ge_ABS
Follow this repository to see the codes related to simulating the spin splitting of Andreev Bound Spectrum with germanium as the weak link.

Following files are important and contain the models
1. Operators2D.py
2. Minimal_model_1D.py
3. Minimal_model_1D_1sub.py
4. Minimal_model_2D.py(has disorder included in it)
5. charge_trap_pot_model.py
6. charge_trap_pot_model_neg_charge.py

All the required libraries and global constants are imported in each of the above scipt using 
5. config.py file

This file is independent of all other file
6. Eff_heavyhole_model_CSOC_parameters.py (is used to find out Ge parameters using eff heavy hole model) 
7. minimal_model_analysis_with_W_mu_subbands_LCSOC_notebooks/Eff_mass_and_SO_vs_Lz.ipynb - The actual parameters used in the paper are found in this file. 

Ignore the below files. Refer to chtc_Spin_Splitting_Ge_ABS repo for correct version of them. Didnt use any of these in the paper
8. spin_splitting_map_cubic.py
9. spin_splitting_map_linear.py  

10. venv is the virtual enviroment
11. Requiements.txt lists all necessary libraries except the libraries required to run jupyter notebooks in VS code.
12. If you need to rum Jupyter notebooks, install the relevant packages for that too, if not found after requirements.txt

These are the analysis files for the manuscript. 
13. Fig1 - Figure_code_final/Andreev_Bound_Spectrum.ipynb
14. Fig2 - Figure_Codes\Variation_ss_soc_linear_cubic_2d_together_75nm.ipynb
15. Fig3 - Figure_code_final\lsoc_csoc_spectrum.ipynb, Figure_code_final\lsoc_csoc_spectrum_zoom.ipynb
16. Fig4 - Figure_code_final\Es_vs_omega_mu__zero_gradient_csoc_lsoc_together_final.ipynb,                 disorder_charge_trap\cubic_linear_spin_Splitting_vs_yct_analysis_together.ipynb 
17. Appendix D - sheet_of_charge_v2.ipynb, Image_charge_3d.py

