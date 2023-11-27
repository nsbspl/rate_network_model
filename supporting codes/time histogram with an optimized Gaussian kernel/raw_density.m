% spk_trains = F_binary_5Hz_Vim;
% raw_density = raw_density(F_binary_5Hz_Vim);

function final_raw_density = raw_density(spk_trains)

% concatenate all the spike trains recorded for a certain DBS fq
% and compute the raw density as the : conc_spk_train/"num of trials"


len = size(spk_trains,2);
num_trial = size(spk_trains,1);
conc_spk_train = zeros(1,len);

for i = 1:num_trial
    spk_this_trial = spk_trains(i,:);
    conc_spk_train = conc_spk_train + spk_this_trial;
end

final_raw_density = conc_spk_train/num_trial;

end

%% Compare the conc spk train with a single spk train
% figure; hold on
% ax(1) = subplot(2,1,1); plot(spk_trains(1,:),'r')
% ax(2) = subplot(2,1,2); plot(final_raw_density)
% linkaxes(ax,'x')
% hold off

% figure; hold on
% ax(1) = subplot(2,1,1); plot(x,'r')
% ax(2) = subplot(2,1,2); plot(x_ab)
% linkaxes(ax,'x')
% hold off

% x = conc_spk_train;


