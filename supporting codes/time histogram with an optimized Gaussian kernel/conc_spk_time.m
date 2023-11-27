% concatenate all the spike timing (in sec) in multiple spike train trials 
% to compute the Hideaki PSTH_Gau kernel


% dt = 0.1; (ms)
% spk_trains = F_binary_100Hz_Vim;

function conc_spk_timing = conc_spk_time(spk_trains, dt)

num_trials = size(spk_trains,1); % # of spike train trials 
idx_spk = []; % collect spike indices

for i = 1:num_trials
    idx_spk_this_trial = find(spk_trains(i,:) == 1);
    idx_spk = [idx_spk,idx_spk_this_trial];
end

idx_spk = sort(idx_spk); %  sorts the spike indices in ascending order.

conc_spk_timing = (idx_spk-1)*dt*1e-3; % change indices into time (s)


end

