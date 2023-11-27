function NMSE = normalized_MSE(signal, pred)
err_vec = signal - pred;
root_MSE = sqrt(sum(err_vec.^2)/length(err_vec));
 %signal = signal - min(signal);
RSS = sqrt(sum(signal.^2)/length(signal));
NRMSE = root_MSE/RSS;
NMSE = NRMSE^2;

% err_vec = err_ex1;
% signal = baseline_ex1_0;
% NRMSE_ex102 = normalized_rootMSE(baseline_ex102_0,baseline_ex102_err);
% NRMSE_fast_vec = [NRMSE_fast_vec, NRMSE_ex97];
% NRMSE_mixed_vec = [NRMSE_mixed_vec, NRMSE_ex101];
% NRMSE_slow_vec = [NRMSE_slow_vec, NRMSE_ex80];

% signal = FR_psth_50Hz_SNr; pred =opt_fit_FR_50Hz_SNr; NMSE = normalized_MSE(signal, pred)



% mean_NRMSE_fast = mean(NRMSE_fast_vec);
% sd_NRMSE_fast = sqrt(var(NRMSE_fast_vec));

%plot(1:length(baseline_ex3_0),baseline_ex3_0,1:length(baseline_ex3_0),lpc_pred_ex3,'--','Linewidth',1); set(gca,'FontSize',20);
%   figure; hold on
%    ax(1) = subplot(3,1,1); plot(baseline_ex3_0)
%    title('fast baseline artifact')
%    ax(2) = subplot(3,1,2); plot(removed_sig_s)
%    title('slow baseline artifact')
%    ax(3) = subplot(3,1,3); plot(removed_sig_m)
%    title('mixed baseline artifact')
%    linkaxes(ax,'x')
%    hold off