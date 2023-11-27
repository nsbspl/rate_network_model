function [opt_fit_FR_f,para_vec, sim_r_e_f, sim_r_i_f]  = Stabilizing_Function_Refining_Stage(ini_para_vec, r_D_0, r_I_0, FR_time_hist,I_DBS)

%% Introduction

% Stabilizing Function in Refining Stage

% FR_time_hist: Experimental firing rate computed by time histogram. 
% Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated.

% I_DBS: the DBS-induced post-synaptic current, obtained by the Tsodyks & Markram model (Tsodyks et al. (1998)). 
% Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated.

% ini_para_vec = [tau_i, Wee, Wie, Wei, Wii, tau_E, c, s, k, r_E_b]; initial parameters
% tau_i: inhibitory time constant; tau_E: excitatory time constant
% Wee, Wie, Wei, Wii: connectivity strengths (see paper)
% c, s, k: sigmoid parameters (see paper)
% r_E_b (variable): baseline firing rate of the external excitatory neural group

% r_I_0 (constant): baseline firing rate of the inhibitory neural group
% r_D_0 (constant): baseline firing rate of Vim neurons


% rate network model consists of 3 ODEs:
% (1). The Vim neural group receiving DBS: r
% (2). The external exc neural group: r_e
% (3). The external inh neural group: r_i

%% inputs

% r_D_0 = 25;  r_I_0 = 5; FR_time_hist = FR_time_hist_conc; I_DBS = I_DBS_conc;

%% 1. Initial and optimization settings
dt = 0.1; %ms
sim_r = zeros(length(FR_time_hist),1);
sim_r_e = zeros(length(FR_time_hist),1);
sim_r_i = zeros(length(FR_time_hist),1);

sim_r(1) = r_D_0;
sim_r_e(1) = ini_para_vec(10);
sim_r_i(1) = r_I_0;

dr = zeros(length(FR_time_hist),1);
dr_e = zeros(length(FR_time_hist),1);
dr_i = zeros(length(FR_time_hist),1);

delay_used = 0; % delay related to time histogram computing (not needed if data appropriately processed)

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

delay = delay_used; 
% Note that this variable is repeated in SSE_fn. If delay_used is non-zero, adjust it within SSE_fn to
% be consistent
delay_idx = floor(delay/dt);


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
    
    Wr_e = para_vec(2)*(sim_r(i-1) + sim_r_e(i-1)) - para_vec(4)*sim_r_i(i-1);
    Wr_i = para_vec(3)*(sim_r(i-1) + sim_r_e(i-1)) - para_vec(5)*sim_r_i(i-1);
    
    dr(i) = (dt/para_vec(6))*(Wr_e-(sim_r(i-1)-r_D_0)+ para_vec(7)/(1+exp(-para_vec(8)*(I_DBS(i)-para_vec(9)))));
    dr_e(i) = (dt/para_vec(6))*(Wr_e-(sim_r_e(i-1)-para_vec(10)));
    dr_i(i) = (dt/para_vec(1))*(Wr_i-(sim_r_i(i-1)-r_I_0));
    
        
    sim_r(i) = sim_r(i-1) + dr(i);
    sim_r_e(i) = sim_r_e(i-1) + dr_e(i);
    sim_r_i(i) = sim_r_i(i-1) + dr_i(i);
        
    
    if sim_r(i) <= 0
        sim_r(i) = 0;
    end
    
    if sim_r_e(i) <= 0
        sim_r_e(i) = 0;
    end
    
    if sim_r_i(i) <= 0
        sim_r_i(i) = 0;
    end
    
    if i == end_idx_5Hz+1 || i == end_idx_10Hz+1 || i == end_idx_20Hz+1 || i == end_idx_30Hz+1 || i == end_idx_50Hz+1 || i == end_idx_100Hz+1 
        sim_r(i) =  r_D_0; 
        sim_r_e(i) = para_vec(10);
        sim_r_i(i) = r_I_0;
        
        % These time points correspond to the start of a new fq DBS, thus the sim_FR is reset to the steady state.
        % This is exactly equivalent to "seperating the signal w.r.t fq segments,
        % and sum the SSE for each segment".
    end
    
end

SSE_now = 0;

delay = delay_used; % delay related to time histogram computing (not needed if data appropriately processed)
delay_idx = floor(delay/dt);





%%  weight of each DBS frequency

  ratio_5 = 2000; ratio_10 = 100;ratio_20 = 400;ratio_30 = 500;ratio_50 = 3000;ratio_100 = 12000;ratio_200 = 150000; 

% match the # of pulses (each DBS fq data has the same # of pulses)
num_5 = 480; num_10 = 480; num_20 = 400;num_30 = 400;num_50 = 480;num_100 = 48;num_200 = 60;

weight_5Hz = ratio_5*num_5;
weight_10Hz = ratio_10*num_10;
weight_20Hz = ratio_20*num_20;
weight_30Hz = ratio_30*num_30;
weight_50Hz = ratio_50*num_50;
weight_100Hz = ratio_100*num_100;
weight_200Hz = ratio_200*num_200;



for i = 1:(length(FR_time_hist)-delay_idx)
    if i < end_idx_5Hz+1-delay_idx
        SSE_now = SSE_now + weight_5Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
    
    if i < end_idx_10Hz+1-delay_idx && i > end_idx_5Hz
        SSE_now = SSE_now + weight_10Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
    
    if i < end_idx_20Hz+1-delay_idx && i > end_idx_10Hz
        SSE_now = SSE_now + weight_20Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
    
    if i < end_idx_30Hz+1-delay_idx && i > end_idx_20Hz
        SSE_now = SSE_now + weight_30Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
    
    if i < end_idx_50Hz+1-delay_idx && i > end_idx_30Hz
        SSE_now = SSE_now + weight_50Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
    
    if i < end_idx_100Hz+1-delay_idx && i > end_idx_50Hz
        SSE_now = SSE_now + weight_100Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
    
    if i < end_idx_200Hz+1-delay_idx && i > end_idx_100Hz
        SSE_now = SSE_now + weight_200Hz*(sim_r(i)-FR_time_hist(i+delay_idx))^2;
    end
end



%% Constraints 

e_constant_upper_limit = 9.4; 

e_constant_lower_limit = 0.5; 

i_constant_upper_limit = 24; 

i_constant_lower_limit = 5; 

baseline_FR_E_lower_limit = 10;

baseline_FR_E_upper_limit = 70; 


if para_vec(1) < 2*para_vec(6) || para_vec(1) < i_constant_lower_limit || para_vec(2) < 0.2 || para_vec(3) < 0 || para_vec(4) < 0 ||  para_vec(5) < 0 || para_vec(6)< e_constant_lower_limit || para_vec(1) > i_constant_upper_limit || para_vec(6)> e_constant_upper_limit  || para_vec(10) < baseline_FR_E_lower_limit || para_vec(10) > baseline_FR_E_upper_limit || para_vec(2) > para_vec(4) || para_vec(3) > para_vec(4) || para_vec(5) > para_vec(4) || para_vec(5) < 0.001 || para_vec(3) < 0.001
    SSE_now = 1e222*SSE_now;
end

%  W_ie and W_ii are both small

if para_vec(5) < 1 
    SSE_now = 1e-2*SSE_now;
end

if para_vec(5) < 0.5 
    SSE_now = 1e-2*SSE_now;
end

if para_vec(3) < 0.05 
    SSE_now = 1e-2*SSE_now;
end

if para_vec(5) < 0.05 
    SSE_now = 1e-2*SSE_now;
end


end

%% 3. Optimal simulation & results

% A. simulate again with the optimal para's

opt_fit_FR = zeros(length(FR_time_hist),1);
sim_r_e = zeros(length(FR_time_hist),1);
sim_r_i = zeros(length(FR_time_hist),1);

opt_fit_FR(1) = r_D_0;
sim_r_e(1) = ini_para_vec(10);
sim_r_i(1) = r_I_0;

dr = zeros(length(FR_time_hist),1);
dr_e = zeros(length(FR_time_hist),1);
dr_i = zeros(length(FR_time_hist),1);



for i = 2:length(FR_time_hist)
    
    Wr_e = para_vec(2)*(opt_fit_FR(i-1) + sim_r_e(i-1)) - para_vec(4)*sim_r_i(i-1);
    Wr_i = para_vec(3)*(opt_fit_FR(i-1) + sim_r_e(i-1)) - para_vec(5)*sim_r_i(i-1);
    
    dr(i) = (dt/para_vec(6))*(Wr_e-(opt_fit_FR(i-1)-r_D_0)+para_vec(7)/(1+exp(-para_vec(8)*(I_DBS(i)-para_vec(9)))));
    dr_e(i) = (dt/para_vec(6))*(Wr_e-(sim_r_e(i-1)-para_vec(10)));
    dr_i(i) = (dt/para_vec(1))*(Wr_i-(sim_r_i(i-1)-r_I_0));
           
    opt_fit_FR(i) = opt_fit_FR(i-1) + dr(i);
    sim_r_e(i) = sim_r_e(i-1) + dr_e(i);
    sim_r_i(i) = sim_r_i(i-1) + dr_i(i);
        
    
    if opt_fit_FR(i) <= 0
        opt_fit_FR(i) = 0;
    end
    
    if sim_r_e(i) <= 0
        sim_r_e(i) = 0;
    end
    
    if sim_r_i(i) <= 0
        sim_r_i(i) = 0;
    end
    
    if i == end_idx_5Hz+1 || i == end_idx_10Hz+1 || i == end_idx_20Hz+1 || i == end_idx_30Hz+1 || i == end_idx_50Hz+1 || i == end_idx_100Hz+1 
        opt_fit_FR(i) =  r_D_0;
        sim_r_e(i) = para_vec(10);
        sim_r_i(i) = r_I_0;
        
        % These time points correspond to the start of a new fq DBS, thus the sim_FR is reset to the steady state.
        % This is exactly equivalent to "seperating the signal w.r.t fq segments,
        % and sum the SSE for each segment".
    end
end

%% Plots

%%  Plots

% match experimental data and model sim (if delay is non-zero)
opt_fit_FR_5Hz = opt_fit_FR(1:end_idx_5Hz-delay_idx);
opt_fit_FR_10Hz = opt_fit_FR(end_idx_5Hz+1:end_idx_10Hz-delay_idx);
opt_fit_FR_20Hz = opt_fit_FR(end_idx_10Hz+1:end_idx_20Hz-delay_idx);
opt_fit_FR_30Hz = opt_fit_FR(end_idx_20Hz+1:end_idx_30Hz-delay_idx);
opt_fit_FR_50Hz = opt_fit_FR(end_idx_30Hz+1:end_idx_50Hz-delay_idx);
opt_fit_FR_100Hz = opt_fit_FR(end_idx_50Hz+1:end_idx_100Hz-delay_idx);
opt_fit_FR_200Hz = opt_fit_FR(end_idx_100Hz+1:end_idx_200Hz-delay_idx);
opt_fit_FR_f = [opt_fit_FR_5Hz;opt_fit_FR_10Hz;opt_fit_FR_20Hz;opt_fit_FR_30Hz;opt_fit_FR_50Hz;opt_fit_FR_100Hz;opt_fit_FR_200Hz];

FR_5Hz = FR_time_hist(1+delay_idx : end_idx_5Hz);
FR_10Hz = FR_time_hist(end_idx_5Hz+1+delay_idx:end_idx_10Hz);
FR_20Hz = FR_time_hist(end_idx_10Hz+1+delay_idx:end_idx_20Hz);
FR_30Hz = FR_time_hist(end_idx_20Hz+1+delay_idx:end_idx_30Hz);
FR_50Hz = FR_time_hist(end_idx_30Hz+1+delay_idx:end_idx_50Hz);
FR_100Hz = FR_time_hist(end_idx_50Hz+1+delay_idx:end_idx_100Hz);
FR_200Hz = FR_time_hist(end_idx_100Hz+1+delay_idx:end_idx_200Hz);
FR_time_hist_f = [FR_5Hz;FR_10Hz;FR_20Hz;FR_30Hz;FR_50Hz;FR_100Hz;FR_200Hz];

% 2. indices 
idx_all = 1:length(opt_fit_FR_f);

idx_5Hz = (1:end_idx_5Hz-delay_idx); 
idx_10Hz = (end_idx_5Hz+1-delay_idx:end_idx_10Hz-2*delay_idx);
idx_20Hz = (end_idx_10Hz+1-2*delay_idx:end_idx_20Hz-3*delay_idx);
idx_30Hz = (end_idx_20Hz+1-3*delay_idx:end_idx_30Hz-4*delay_idx);
idx_50Hz = (end_idx_30Hz+1-4*delay_idx:end_idx_50Hz-5*delay_idx);
idx_100Hz = (end_idx_50Hz+1-5*delay_idx:end_idx_100Hz-6*delay_idx);
idx_200Hz = (end_idx_100Hz+1-6*delay_idx:end_idx_200Hz-7*delay_idx);

 % (8). All together
figure(8)
plot((idx_all-1)*dt,FR_time_hist_f,'g', (idx_all-1)*dt,opt_fit_FR_f,'--','Linewidth',1); 
set(gca,'FontSize',20)
title('Experimental firing rate and network model simulation, Vim, all DBS frequencies (5~200)')
legend({'experimental firing rate','network model simulation'},'FontSize', 24)

% sim_r_e and sim_r_i

sim_r_e_5Hz = sim_r_e(1:end_idx_5Hz-delay_idx);
sim_r_e_10Hz = sim_r_e(end_idx_5Hz+1:end_idx_10Hz-delay_idx);
sim_r_e_20Hz = sim_r_e(end_idx_10Hz+1:end_idx_20Hz-delay_idx);
sim_r_e_30Hz = sim_r_e(end_idx_20Hz+1:end_idx_30Hz-delay_idx);
sim_r_e_50Hz = sim_r_e(end_idx_30Hz+1:end_idx_50Hz-delay_idx);
sim_r_e_100Hz = sim_r_e(end_idx_50Hz+1:end_idx_100Hz-delay_idx);
sim_r_e_200Hz = sim_r_e(end_idx_100Hz+1:end_idx_200Hz-delay_idx);
sim_r_e_f = [sim_r_e_5Hz;sim_r_e_10Hz;sim_r_e_20Hz;sim_r_e_30Hz;sim_r_e_50Hz;sim_r_e_100Hz;sim_r_e_200Hz];

sim_r_i_5Hz = sim_r_i(1:end_idx_5Hz-delay_idx);
sim_r_i_10Hz = sim_r_i(end_idx_5Hz+1:end_idx_10Hz-delay_idx);
sim_r_i_20Hz = sim_r_i(end_idx_10Hz+1:end_idx_20Hz-delay_idx);
sim_r_i_30Hz = sim_r_i(end_idx_20Hz+1:end_idx_30Hz-delay_idx);
sim_r_i_50Hz = sim_r_i(end_idx_30Hz+1:end_idx_50Hz-delay_idx);
sim_r_i_100Hz = sim_r_i(end_idx_50Hz+1:end_idx_100Hz-delay_idx);
sim_r_i_200Hz = sim_r_i(end_idx_100Hz+1:end_idx_200Hz-delay_idx);
sim_r_i_f = [sim_r_i_5Hz;sim_r_i_10Hz;sim_r_i_20Hz;sim_r_i_30Hz;sim_r_i_50Hz;sim_r_i_100Hz;sim_r_i_200Hz];

% 

% 9. plot with spk time
% figure(8)
% plot(sim_r_i)
% plot(sim_r_e)

end

