% plot_data

mean_abs_pct_diff_ini_vec = [mean_abs_pct_diff_ini_2;mean_abs_pct_diff_ini_4;mean_abs_pct_diff_ini_6;mean_abs_pct_diff_ini_10;mean_abs_pct_diff_ini_14;mean_abs_pct_diff_ini_18;mean_abs_pct_diff_ini_22;mean_abs_pct_diff_ini_26;mean_abs_pct_diff_ini_30;mean_abs_pct_diff_ini_34;mean_abs_pct_diff_ini_38;mean_abs_pct_diff_ini_42;mean_abs_pct_diff_ini_46;mean_abs_pct_diff_ini_50];

mean_abs_pct_diff_final_vec = [mean_abs_pct_diff_final_2;mean_abs_pct_diff_final_4;mean_abs_pct_diff_final_6;mean_abs_pct_diff_final_10;mean_abs_pct_diff_final_14;mean_abs_pct_diff_final_18;mean_abs_pct_diff_final_22;mean_abs_pct_diff_final_26;mean_abs_pct_diff_final_30;mean_abs_pct_diff_final_34;mean_abs_pct_diff_final_38;mean_abs_pct_diff_final_42;mean_abs_pct_diff_final_46;mean_abs_pct_diff_final_50];

ER_diff_2 = ER_2 - ER_ref;
ER_diff_4 = ER_4 - ER_ref;ER_diff_6 = ER_6 - ER_ref;ER_diff_10 = ER_10 - ER_ref;ER_diff_14 = ER_14 - ER_ref;
ER_diff_18 = ER_18 - ER_ref;ER_diff_22 = ER_22 - ER_ref;ER_diff_26 = ER_26 - ER_ref;ER_diff_30 = ER_30 - ER_ref;
ER_diff_34 = ER_34 - ER_ref;ER_diff_38 = ER_38 - ER_ref;ER_diff_42 = ER_42 - ER_ref;ER_diff_46 = ER_46 - ER_ref;
ER_diff_50 = ER_50 - ER_ref;

ER_vec = [ER_2;ER_4;ER_6;ER_10;ER_14;ER_18;ER_22;ER_26;ER_30;ER_34;ER_38;ER_42;ER_46;ER_50];
rou_D_vec = [rou_D_2;rou_D_4;rou_D_6;rou_D_10;rou_D_14;rou_D_18;rou_D_22;rou_D_26;rou_D_30;rou_D_34;rou_D_38;rou_D_42;rou_D_46;rou_D_50];
ER_diff_vec = [ER_diff_2;ER_diff_4;ER_diff_6;ER_diff_10;ER_diff_14;ER_diff_18;ER_diff_22;ER_diff_26;ER_diff_30;ER_diff_34;ER_diff_38;ER_diff_42;ER_diff_46;ER_diff_50];




% Example vectors (replace with your own)
x = 100*mean_abs_pct_diff_ini_vec;
%y = 100*mean_abs_pct_diff_final_vec;
%y = 100*ER_vec;
%y = 100*ER_abs_diff_vec;
y = rou_D_vec;


% Color gradient variable (can be x, y, or index)

%c = linspace(1,length(x),length(x));
c = [2,4,6,10,14,18,22,26,30,34,38,42,46,50];

figure;

scatter(x, y, 60, c, 'filled');  % 60 = marker size

colormap('jet');    % choose colormap: parula, jet, turbo, viridis
cb = colorbar;            % optional color scale
cb.Ticks = [2,6,10,14,18,22,26,30,34,38,42,46,50];
%set(gca,'FontSize',20,'LineWidth',1);

xlabel('|diff| in initial parameters (%)');
%ylabel('|diff| in final parameters (%)');
ylabel('rou_D');
%title('compare mean absolute percentage difference in parameters');
ylim([0,1])

set(gca,'FontSize',20,'LineWidth',1);