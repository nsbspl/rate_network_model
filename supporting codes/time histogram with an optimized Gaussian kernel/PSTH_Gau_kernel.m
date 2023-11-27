% compute the PSTH with a Gaussian kernel, e.g., the optimized time
% histogram Gaussian kernel from the Hideaki method

% inputs:
% (1).SPK_all_trial: all spike trains, 
% # of rows = # of spike train trials; # of cols = # of sampling points
% (2). Gau_kernel: the time histogram Gaussian kernel in sec 
% (3). dt = 0.1; %ms, sampling resolution.


% SPK_all_trial = F_binary_100Hz_Vim; Gau_kernel = optw; dt = 0.1;


function PSTH_FR_Gau = PSTH_Gau_kernel(SPK_all_trial, Gau_kernel,dt)

TW = 1000*Gau_kernel; %ms

% psth_total = PSTH_my(NG)(SPK_all_trial', dt, dt); % same as "conc_spk_trains" in the next line 
conc_spk_trains = sum(SPK_all_trial,1)'; % "the concatenated spike trains of all trials"

Trial_num = size(SPK_all_trial,1);

FG_mix = KernelPSTH_ori(conc_spk_trains,TW,dt);
fr_neuron = sum(conc_spk_trains)/((size(SPK_all_trial,2)-1)*dt*1e-3)/Trial_num;
scale_ = fr_neuron/mean(FG_mix); % sum(scale_*FG_mix)*p.dt*1e-3/EndTime = sum(psth_total) / Trial_num
FG_mix = scale_*FG_mix;

PSTH_FR_Gau = FG_mix;

plot(PSTH_FR_Gau)

end