%

end_idx_5Hz = 96362;
end_idx_10Hz = 143958;
end_idx_20Hz = 172542;
end_idx_30Hz = 191876;
end_idx_50Hz = 201077;
end_idx_100Hz = 251078;
end_idx_200Hz = 270779;

test_sig = final_fit;
%test_sig = attempt_fit_FR;


NMSE_5_to_50Hz = normalized_MSE(PSTH_FR_5_to_50_Hz,test_sig(1:end_idx_50Hz));
%% NMSE
NMSE_5Hz = normalized_MSE(FR_PSTH_5Hz,test_sig(1:end_idx_5Hz));
NMSE_10Hz = normalized_MSE(FR_PSTH_10Hz,test_sig(end_idx_5Hz+1:end_idx_10Hz));
NMSE_20Hz = normalized_MSE(FR_PSTH_20Hz,test_sig(end_idx_10Hz+1:end_idx_20Hz));
NMSE_30Hz = normalized_MSE(FR_PSTH_30Hz,test_sig(end_idx_20Hz+1:end_idx_30Hz));
NMSE_50Hz = normalized_MSE(FR_PSTH_50Hz,test_sig(end_idx_30Hz+1:end_idx_50Hz));

NMSE_100Hz = normalized_MSE(FR_PSTH_100Hz,test_sig(end_idx_50Hz+1:end_idx_100Hz));

NMSE_200Hz = normalized_MSE(FR_PSTH_200Hz,test_sig(end_idx_100Hz+1:end_idx_200Hz));


NMSE_vec = [NMSE_5_to_50Hz,NMSE_100Hz,NMSE_200Hz];
ER = mean(NMSE_vec);

%NMSE_vec = [NMSE_5Hz,NMSE_10Hz,NMSE_20Hz,NMSE_30Hz,NMSE_50Hz];

% 
% mean(NMSE_vec)

%% R^2
% NMSE_5Hz = normalized_MSE(FR_PSTH_5Hz,test_sig(1:end_idx_5Hz));
% NMSE_10Hz = normalized_MSE(FR_PSTH_10Hz,test_sig(end_idx_5Hz+1:end_idx_10Hz));
% NMSE_20Hz = normalized_MSE(FR_PSTH_20Hz,test_sig(end_idx_10Hz+1:end_idx_20Hz));
% NMSE_30Hz = normalized_MSE(FR_PSTH_30Hz,test_sig(end_idx_20Hz+1:end_idx_30Hz));
% NMSE_50Hz = normalized_MSE(FR_PSTH_50Hz,test_sig(end_idx_30Hz+1:end_idx_50Hz));
% 
% NMSE_100Hz = normalized_MSE(FR_PSTH_100Hz,test_sig(end_idx_50Hz+1:end_idx_100Hz));
% 
% NMSE_200Hz = normalized_MSE(FR_PSTH_200Hz,test_sig(end_idx_100Hz+1:end_idx_200Hz));
