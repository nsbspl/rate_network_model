dt = 0.1; %ms


% idx_5Hz = 1:length(FR_5Hz);
% idx_10Hz = 1:length(FR_10Hz);
% idx_20Hz = 1:length(FR_20Hz);
% idx_30Hz = 1:length(FR_30Hz);
% idx_50Hz = 1:length(FR_50Hz);
% idx_100Hz = 1:length(FR_100Hz);
% idx_200Hz = 1:length(FR_200Hz);

% end_idx_5Hz = 96362;
% end_idx_10Hz = 143958;
% end_idx_20Hz = 172542;
% end_idx_30Hz = 191876;
% end_idx_50Hz = 201077;
% end_idx_100Hz = 251078;
% end_idx_200Hz = 270779;
% 
% delay_idx = 0;
% 
% fit_FR_5Hz = final_fit(1+delay_idx : end_idx_5Hz);
% fit_FR_10Hz = final_fit(end_idx_5Hz+1+delay_idx:end_idx_10Hz);
% fit_FR_20Hz = final_fit(end_idx_10Hz+1+delay_idx:end_idx_20Hz);
% fit_FR_30Hz = final_fit(end_idx_20Hz+1+delay_idx:end_idx_30Hz);
% fit_FR_50Hz = final_fit(end_idx_30Hz+1+delay_idx:end_idx_50Hz);
% fit_FR_100Hz = final_fit(end_idx_50Hz+1+delay_idx:end_idx_100Hz);
% fit_FR_200Hz = final_fit(end_idx_100Hz+1+delay_idx:end_idx_200Hz);


ref_FR = FR_200Hz;
fit_FR = fit_FR_200Hz;
idx_used = idx_200Hz;


figure(8)
plot((idx_used-1)*dt,ref_FR,'k', (idx_used-1)*dt,fit_FR,'g','Linewidth',2);
set(gca,'FontSize',24,'Linewidth',1)
ylim([0,500])
xlim([0,2000])

title('clinical firing rate and network model simulation, Vim')
legend({'clinical firing rate','network model simulation'},'FontSize', 24)