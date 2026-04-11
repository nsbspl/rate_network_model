
% ini_para_vec = [tau_i, Wee, Wie, Wei, Wii, tau_E, c, s, k, r_E_b]; 

% NMSE Data


values = 100*[NMSE_5Hz NMSE_10Hz NMSE_20Hz NMSE_30Hz NMSE_50Hz NMSE_100Hz NMSE_200Hz];

% Category labels
categories = {'5Hz','10Hz','20Hz','30Hz','50Hz','100Hz','200Hz'};

% Create bar plot
figure;
b = bar(values);

% Adjust bar width
b.BarWidth = 0.5; 

% Set x-axis labels
set(gca,'XTick',1:7,'XTickLabel',categories);

% Axis labels
xlabel('DBS frequencies');
ylabel('NMSE(%)');
%yticks([0 10 20 30 40 50])

% Optional formatting
set(gca,'FontSize',24,'LineWidth',1);

% Grid
% grid on;