%% Introduction
% route optimization in Global Stage

% FR_time_hist: Experimental firing rate computed by time histogram. 
% Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated.

% I_DBS: the DBS-induced post-synaptic current, obtained by the Tsodyks & Markram model (Tsodyks et al. (1998)). 
% Data from all DBS frequencies {5,10,20,30,50,100,200Hz} are concatenated.

% ini_para_vec: initial parameters going into Global Stage

%% inputs

FR_time_hist = FR_time_hist_conc; I_DBS = I_DBS_conc;
r_I_0 = 5; % baseline firing rate of the inhibitory neural group
r_D_0 = 25; % baseline firing rate of Vim neurons
num_loop = 4; % # of loops in the sequential optimization

initial_para = para_vec; % "para_vec" is obtained from the preliminary fit (see Readme of this folder)

%% optimization

delay = 0; % delay related to time histogram computing (not needed if data appropriately processed)
dt = 0.1; %ms
delay_idx = floor(delay/dt);
sim_length = length(FR_time_hist) - 7*delay_idx; % actual sim data length after removing delays


opt_fit_mx = zeros(sim_length, 2*num_loop);
NMSE_vec = zeros(1,2*num_loop);
para_mx = zeros(length(initial_para),2*num_loop);


para_vec_S = initial_para;
for i = 1:num_loop
    
    [opt_fit_W,para_vec_P, ~, ~]  = Pushing_Function_Global_Stage(para_vec_S, r_D_0, r_I_0,FR_time_hist,I_DBS);
    opt_fit_mx(:,2*i-1) = opt_fit_W;
    NMSE_vec(2*i-1) = normalized_MSE(FR_time_hist, opt_fit_W);
    para_mx(:,2*i-1) = para_vec_P;
    % para_vec_P is the result obtained from a Pushing Function

    
    [opt_fit_AP,para_vec_S, sim_r_e, sim_r_i]  = Stabilizing_Function_Global_Stage(para_vec_P, r_D_0, r_I_0,FR_time_hist,I_DBS);
    opt_fit_mx(:,2*i) = opt_fit_AP;
    NMSE_vec(2*i) = normalized_MSE(FR_time_hist, opt_fit_AP);
    para_mx(:,2*i) = para_vec_S;
    % para_vec_S is the result obtained from a Stabilizing Function
    

end


%% Plots
f_fit = opt_fit_mx(:,end);
idx_all = 1:length(f_fit);


end_idx_5Hz = 96362;
end_idx_10Hz = 143958;
end_idx_20Hz = 172542;
end_idx_30Hz = 191876;
end_idx_50Hz = 201077;
end_idx_100Hz = 251078;
end_idx_200Hz = 270779;

FR_5Hz = FR_time_hist(1+delay_idx : end_idx_5Hz);
FR_10Hz = FR_time_hist(end_idx_5Hz+1+delay_idx:end_idx_10Hz);
FR_20Hz = FR_time_hist(end_idx_10Hz+1+delay_idx:end_idx_20Hz);
FR_30Hz = FR_time_hist(end_idx_20Hz+1+delay_idx:end_idx_30Hz);
FR_50Hz = FR_time_hist(end_idx_30Hz+1+delay_idx:end_idx_50Hz);
FR_100Hz = FR_time_hist(end_idx_50Hz+1+delay_idx:end_idx_100Hz);
FR_200Hz = FR_time_hist(end_idx_100Hz+1+delay_idx:end_idx_200Hz);
FR_time_hist_f = [FR_5Hz;FR_10Hz;FR_20Hz;FR_30Hz;FR_50Hz;FR_100Hz;FR_200Hz];

% 
 % (8). All together
figure(8)
plot((idx_all-1)*dt,FR_time_hist_f,'g', (idx_all-1)*dt,opt_fit_mx(:,end),'--','Linewidth',1);
set(gca,'FontSize',20)
title('Experimental firing rate and network model simulation, Vim, all DBS frequencies (5~200)')
legend({'experimental firing rate','network model simulation'},'FontSize', 24)



% plot(sim_r_i)
% plot(sim_r_e)
    
    
