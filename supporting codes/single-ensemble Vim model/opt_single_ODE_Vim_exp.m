function [opt_fit_FR,para_vec]  = opt_single_ODE_Vim_exp(ini_para_vec,r_0, FR_time_hist,I_DBS)

%% Introduction

% The single-ensemble Vim model in Tian et al. (2023) (Neuromodulation)

% Minimizing SSE to fit the experimental Vim firing rate 
% (computed by a time histogram method with multiple experimental spike trains).

% r_0 = 25 is the baseline firing rate of Vim neurons (DBS-OFF)

% tau_0, c_0, s_0, k_0, h_0, are initial parameters.
% tau_0 is the time constant. c_0, s_0 and k_0 are parameters in the sigmoid function (Tian et al. (2023))
% h is related to r_b in Tian et al. (2023); r_b  = r_0 + h

% FR_time_hist: Experimental firing rate computed by time histogram. 
% Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated.

% I_DBS: the DBS-induced post-synaptic current, obtained by the Tsodyks & Markram model (Tsodyks et al. (1998)). 
% Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated.

%% input parameters

% r_0 = 25; tau_0 = 2; c_0 = 970; s_0 = 4.7e-3; k_0=1900; h_0 = 0;
% ini_para_vec = [tau_0, c_0, s_0, k_0, h_0]; FR_time_hist = FR_time_hist_conc; I_DBS = I_DBS_conc;


%% 1. Initial and optimization settings
dt = 0.1; %ms
sim_FR = zeros(length(FR_time_hist),1);
sim_FR(1) = r_0;
dr = zeros(length(FR_time_hist),1);
para_vec = ini_para_vec;

options = optimset('MaxFunEvals',2e3,'MaxIter',2e3); %initial:'MaxFunEvals'='MaxIter'=1000
para_vec = fminsearch(@SSE_fn,para_vec,options);

end_idx_5Hz = 96362;
end_idx_10Hz = 143958;
end_idx_20Hz = 172542;
end_idx_30Hz = 191876;
end_idx_50Hz = 201077;
end_idx_100Hz = 251078;
end_idx_200Hz = 270779;

%% 2. The objective fn
function SSE_now = SSE_fn(para_vec)

end_idx_5Hz = 96362;
end_idx_10Hz = 143958;
end_idx_20Hz = 172542;
end_idx_30Hz = 191876;
end_idx_50Hz = 201077;
end_idx_100Hz = 251078;
end_idx_200Hz = 270779;

for i = 2:length(FR_time_hist)
    dr(i) = (dt/para_vec(1))*(-(sim_FR(i-1)-r_0-para_vec(5))+para_vec(2)/(1+exp(-para_vec(3)*(I_DBS(i)-para_vec(4)))));
    sim_FR(i) = sim_FR(i-1) + dr(i);
    if sim_FR(i) <= 0
        sim_FR(i) = 0;
    end
    if i == end_idx_5Hz+1 || i == end_idx_10Hz+1 || i == end_idx_20Hz+1 || i == end_idx_30Hz+1 || i == end_idx_50Hz+1 || i == end_idx_100Hz+1 
        sim_FR(i) =  r_0; 
        % These time points correspond to the start of a new fq DBS, thus the sim_FR is reset to the steady state.
        % This is exactly equivalent to "seperating the signal w.r.t fq segments,
        % and sum the SSE for each segment".
    end
    
end
SSE_now = 0;

delay = 0; % delay related to time histogram computing (not needed if data appropriately processed)
delay_idx = floor(delay/dt);


% weight of data of each DBS frequency 
% Note that the recording data length from each DBS frequency is also different.
weight_5Hz = 1;
weight_10Hz = 1;
weight_20Hz = 1;
weight_30Hz = 1;
weight_50Hz = 1;
weight_100Hz = 1;
weight_200Hz = 8;


for i = 1:(length(FR_time_hist)-delay_idx) 
    if i < end_idx_5Hz+1
        SSE_now = SSE_now + weight_5Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
    
    if i < end_idx_10Hz+1 && i > end_idx_5Hz
        SSE_now = SSE_now + weight_10Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
    
    if i < end_idx_20Hz+1 && i > end_idx_10Hz
        SSE_now = SSE_now + weight_20Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
    
    if i < end_idx_30Hz+1 && i > end_idx_20Hz
        SSE_now = SSE_now + weight_30Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
    
    if i < end_idx_50Hz+1 && i > end_idx_30Hz
        SSE_now = SSE_now + weight_50Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
    
    if i < end_idx_100Hz+1 && i > end_idx_50Hz
        SSE_now = SSE_now + weight_100Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
    
    if i < end_idx_200Hz+1 && i > end_idx_100Hz
        SSE_now = SSE_now + weight_200Hz*(sim_FR(i)-FR_time_hist(i+delay_idx))^2; 
    end
end


% exc time constants and the baseline firng rate (r_b = r_0 + h) are restricted within some biological constraints. 

e_constant_upper_limit = 5; 
e_constant_lower_limit = 0.1;
baseline_FR_Vim_upper_limit = 100;
baseline_FR_Vim_lower_limit = 10;

% constraints
if para_vec(1) < e_constant_lower_limit || para_vec(1) > e_constant_upper_limit || para_vec(3) > 0.1 || r_0 + para_vec(5) < baseline_FR_Vim_lower_limit || r_0 + para_vec(5) > baseline_FR_Vim_upper_limit
    SSE_now = 1e222*SSE_now;
end

end

%% 3. Optimal simulation & results

% A. simulate again with the optimal para's
opt_fit_FR = zeros(length(FR_time_hist),1);
opt_fit_FR(1) = r_0;
dr = zeros(length(FR_time_hist),1);
for i = 2:length(FR_time_hist)
    dr(i) = (dt/para_vec(1))*(-(opt_fit_FR(i-1)-r_0-para_vec(5))+para_vec(2)/(1+exp(-para_vec(3)*(I_DBS(i)-para_vec(4)))));
    opt_fit_FR(i) = opt_fit_FR(i-1) + dr(i);
    if opt_fit_FR(i) <= 0
        opt_fit_FR(i) = 0;
    end
    if i == end_idx_5Hz+1 || i == end_idx_10Hz+1 || i == end_idx_20Hz+1 || i == end_idx_30Hz+1 || i == end_idx_50Hz+1 || i == end_idx_100Hz+1 
        opt_fit_FR(i) =  r_0;
        % These time points correspond to the start of a new fq DBS, thus the sim_FR is reset to the steady state.
        % This is exactly equivalent to "seperating the signal w.r.t fq segments,
        % and sum the SSE for each segment".
    end
end

% B. Plots:

idx_all = 1:length(FR_time_hist);
idx_5Hz = (1:end_idx_5Hz); idx_10Hz = (end_idx_5Hz+1:end_idx_10Hz);idx_20Hz = (end_idx_10Hz+1:end_idx_20Hz);idx_30Hz = (end_idx_20Hz+1:end_idx_30Hz);
idx_50Hz = (end_idx_30Hz+1:end_idx_50Hz);idx_100Hz = (end_idx_50Hz+1:end_idx_100Hz);idx_200Hz = (end_idx_100Hz+1:end_idx_200Hz);

 % All together
figure(7)
plot(idx_all*dt,FR_time_hist,'g', idx_all*dt,opt_fit_FR,'--','Linewidth',1); 
set(gca,'FontSize',20)
title('Experimental firing rate and single-ODE model simulation, Vim, all DBS frequencies (5~200)')
legend({'experimental firing rate','single-ODE model simulation'},'FontSize', 24)


end

