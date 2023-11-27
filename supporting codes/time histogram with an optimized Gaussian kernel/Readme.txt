This folder contains the codes for computing an optimized Gaussian time histogram kernel used for computing
instantaneous firing rate from multiple spike trains (Shimazaki et al. (2007) and Shimazaki et al. (2010)). 
Users can try the codes with the provided example F_binary_100Hz_Vim.mat, which contains the experimental spike trains recorded during 100Hz Vim-DBS. Users can implement the following steps:

1. Run conc_spk_time.m to concatenate the spike trains into one vector - conc_spk_timing.mat.
2. Run the code "[y,t,optw] = sskernel_use(conc_spk_timing)" to obtain the optimized kernel in seconds (optw.mat).
3. Run PSTH_Gau_kernel_Y.m to compute the firing rate from the spike trains and the optimized time histogram kernel.