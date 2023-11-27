% Inputs: T = "total simulation time", ST = "starting time of DBS", both in sec.
%         Fq_DBS = "the stimulation frequency of Vim-DBS"

% T=2; ST = 0; Fq_DBS=100; 

% generate the DBS-induced post-synaptic current (I_DBS) using the Tsodyks & Markram model

function I_DBS_total = TM_I_DBS_Vim(T,ST,Fq_DBS)
dt = 0.1; % sampling time (msec)
tt = 0:dt:T*1e3;
L = length(tt);
param.tStop = T*1e3; %(msec)
param.dt = dt;
V = zeros(length(tt),1000);
%% Access to Pre-synaptic neurons (average firing = 7 Hz). Spike times of 1000 neurons (more accurately, 1000 synapses) were saved in Noise_FR.txt 
A_ = dlmread('Noise_FR.txt'); % the rows show the spike time (as it is a matrix, some times are zeros).
%--- Build the spike binary
spTime_indx = A_'/dt;
F_binary = zeros(size(V));
for k=1:size(F_binary,2)
    indx = [];
    indx = find(spTime_indx(:,k)>0 & spTime_indx(:,k)<L);
    F_binary(floor(spTime_indx(indx,k)),k) = 1;
end

%% Setting DBS and Synaptic parameters

ST = ST;
Fs_DBS = Fq_DBS;
activation_percentage = 1; % rate of activated synapses for a given DBS pulse
indx_pool = 1:size(F_binary,2); % 
A = 1;
%--- For excitatory Pre-syanptic neurons
perc_facilitation = 0.4; % percentage of facilitatory synapses % default = 0.5
perc_psudue = 0.2; % percentage of Pseudu Linear synapses % default = 0.2
% 1 - (perc_facilitation + perc_psudue) will be the percentage of
% Depresseion Synapses
Trial_num_Exc = 450; % number of Excitatotry pre-syanptic inputs (out of 500)
Trial_num = Trial_num_Exc;
indx_Exc =indx_pool(randperm(length(indx_pool),Trial_num_Exc)); % randomly selecting indecies from the pool of data
V_sp = F_binary(:,indx_Exc);
tau_syn = 2; % default = 5, opt 22 Mar = 2
Delay = floor((0 + 2*rand(1,Trial_num))/dt); Delay(Delay<0) = 0; % index
tau_f_Fac = 670; tau_d_Fac = 138; U_Fac = 0.19; % For Excitatory & Facilitatory Synapses % default: tau_f_Fac = 670; tau_d_Fac = 138; U_Fac = 0.09;
tau_f_Dep = 17; tau_d_Dep = 85; U_Dep = 0.04; % For Excitatory & Depressing Synapses % default: tau_f_Dep = 17; tau_d_Dep = 671; U_Dep = 0.5;
tau_f_Psu = 326; tau_d_Psu = 329; U_Psu = 0.45;  % For Excitatory & Psudue Linear Synapses % default: tau_f_Psu = 326; tau_d_Psu = 329; U_Psu = 0.29;

run DBS_Plasticity % creates total excitatory current based on parameters of short-term synaptic plasticity
I_Exc = sum(I_cont_delay,2);
%--- For inhibitory Neurons
perc_facilitation = 0.4;
perc_psudue = 0.3;
Trial_num_Inh = 50;
indx_remained = setdiff(indx_pool,indx_Exc);
indx_Inh =indx_remained(randperm(length(indx_remained),Trial_num_Inh));
V_sp = F_binary(:,indx_Inh);
tau_syn = 8.5;% 
Trial_num = Trial_num_Inh;
Delay = floor((0 + 2*rand(1,Trial_num))/dt);  Delay(Delay<0) = 0;% index
tau_f_Fac = 376; tau_d_Fac = 45; U_Fac = 0.016; % For Inhibitory & Facilitatory Synapses
tau_f_Dep = 21; tau_d_Dep = 706; U_Dep = 0.25; % For Inhibitory & Depressing Synapses % default tau_f_Dep = 21; tau_d_Dep = 706; U_Dep = 0.25;
tau_f_Psu = 62; tau_d_Psu = 144; U_Psu = 0.29; % For Inhibitory & Psudue Linear Synapses
run DBS_Plasticity

I_Inh = sum(I_cont_delay,2);%sum(I_cont,2);
w_exc = 3.5; % default: 2.5;
w_inh = 6.0 ; % default: 6.0;
I_DBS_total = 15.5*(w_exc*I_Exc - w_inh*I_Inh); % Total synaptic current % default: 15.0;


plot(I_DBS_total)
