The input data include:

1. FR_time_hist_conc.mat: The experimental firing rate computed by the time histogram method in Shimazaki et al. (2007) and Shimazaki et al. (2010). Data from all Vim-DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated. See the folder "supporting codes/time histogram with an optimized Gaussian kernel" for the codes and an example. 

2. I_DBS_conc.mat: The DBS-induced post-synaptic current, obtained by the Tsodyks & Markram model (Tsodyks et al. (1998) STAR Methods S3 in this paper). Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated. See the folder "supporting codes/Tsodyks & Markram model of STP" for the codes to replicate the results in this paper.

3. ini_para_network_model.mat: The initial parameters, obtained from anatomical literature (see STAR Methods S4.3 in paper for details) and a one-step simulation from our previous single-ensemble model (Tian et al. (2023)). 
Users can run opt_single_ODE_Vim_exp.mat in the folder "supporting codes/single-ensemble Vim model" to obtain the
initial sigmoid parameters {c,s,k}.