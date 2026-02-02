# Spin_Splitting_Ge_ABS <br>
Follow this repository to see the codes related to simulating the spin splitting of Andreev Bound Spectrum with germanium as the weak link. <br>

Following files are important and contain the models <br>
1. Operators2D.py <br>
2. Minimal_model_1D.py <br>
3. Minimal_model_1D_1sub.py <br>
4. Minimal_model_2D.py(has disorder included in it)<br>
5. charge_trap_pot_model.py<br>
6. charge_trap_pot_model_neg_charge.py<br>

All the required libraries and global constants are imported in each of the above scipt using  <br>
5. config.py file <br>

This file is independent of all other file <br>
6. Eff_heavyhole_model_CSOC_parameters.py (is used to find out Ge parameters using eff heavy hole model)  <br>
7. minimal_model_analysis_with_W_mu_subbands_LCSOC_notebooks/Eff_mass_and_SO_vs_Lz.ipynb - The actual parameters used in the paper are found in this file.  <br>

Ignore the below files. Refer to chtc_Spin_Splitting_Ge_ABS repo for correct version of them. Didnt use any of these in the paper <br>
8. spin_splitting_map_cubic.py <br>
9. spin_splitting_map_linear.py   <br>

10. venv is the virtual enviroment <br>
11. Requiements.txt lists all necessary libraries except the libraries required to run jupyter notebooks in VS code. <br>
12. If you need to rum Jupyter notebooks, install the relevant packages for that too, if not found after requirements.txt <br>

These are the analysis files for the manuscript.  <br>
13. Fig1 - Figure_code_final/Andreev_Bound_Spectrum.ipynb <br>
14. Fig2 - Figure_Codes\Variation_ss_soc_linear_cubic_2d_together_75nm.ipynb <br>
15. Fig3 - Figure_code_final\lsoc_csoc_spectrum.ipynb, Figure_code_final\lsoc_csoc_spectrum_zoom.ipynb<br>
16. Fig4 - Figure_code_final\Es_vs_omega_mu__zero_gradient_csoc_lsoc_together_final.ipynb, <br>                                                          disorder_charge_trap\cubic_linear_spin_Splitting_vs_yct_analysis_together.ipynb <br>
17. Appendix D - sheet_of_charge_v2.ipynb, Image_charge_3d.py <br>

