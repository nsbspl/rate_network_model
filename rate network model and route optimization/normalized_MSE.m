function NMSE = normalized_MSE(signal, pred)
err_vec = signal - pred;
root_MSE = sqrt(sum(err_vec.^2)/length(err_vec));
 %signal = signal - min(signal);
RSS = sqrt(sum(signal.^2)/length(signal));
NRMSE = root_MSE/RSS;
NMSE = NRMSE^2;

