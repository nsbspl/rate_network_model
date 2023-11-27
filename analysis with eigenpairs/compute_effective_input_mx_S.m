% compute_effective_input_mx_S

para_vec = opt_para;
fit_r_D = opt_fit_r_D;

mean_r_I = mean(fit_r_i);
mean_r_E = mean(fit_r_e);
mean_r_D = mean(fit_r_D);

s_D_I = mean_r_I * para_vec(4);
s_D_E = (mean_r_E + mean_r_D) * para_vec(2);
s_I_I = mean_r_I * para_vec(5);
s_I_E = (mean_r_E + mean_r_D) * para_vec(3);

S = [s_D_E,s_D_I;s_I_E,s_I_I]

s_D_I_vec = [S_after_loop1(1,2),S_after_loop2(1,2),S_after_loop3(1,2),S_after_loop4(1,2),S_after_loop5(1,2),S_after_loop6(1,2),S_after_loop21(1,2),S_after_loop42(1,2),S_opt_fit(1,2)];
s_D_E_vec = [S_after_loop1(1,1),S_after_loop2(1,1),S_after_loop3(1,1),S_after_loop4(1,1),S_after_loop5(1,1),S_after_loop6(1,1),S_after_loop21(1,1),S_after_loop42(1,1),S_opt_fit(1,1)];

s_I_I_vec = [S_after_loop1(2,2),S_after_loop2(2,2),S_after_loop3(2,2),S_after_loop4(2,2),S_after_loop5(2,2),S_after_loop6(2,2),S_after_loop21(2,2),S_after_loop42(2,2),S_opt_fit(2,2)];
s_I_E_vec = [S_after_loop1(2,1),S_after_loop2(2,1),S_after_loop3(2,1),S_after_loop4(2,1),S_after_loop5(2,1),S_after_loop6(2,1),S_after_loop21(2,1),S_after_loop42(2,1),S_opt_fit(2,1)];


ER_vec = [mean(NMSE_vec_after_loop1),mean(NMSE_vec_after_loop2),mean(NMSE_vec_after_loop3),mean(NMSE_vec_after_loop4),mean(NMSE_vec_after_loop5),mean(NMSE_vec_after_loop6),mean(NMSE_vec_after_loop21),mean(NMSE_vec_after_loop42),mean(NMSE_vec_opt)];

plot(s_D_I_vec,ER_vec,'o')
plot(s_D_E_vec,ER_vec,'o')

plot(s_I_I_vec,ER_vec,'o')
plot(s_I_E_vec,ER_vec,'o')

plot3(s_D_E_vec,s_D_I_vec, ER_vec,'o')

plot3(s_D_I_vec,s_D_E_vec, ER_vec,'o')


% rou_D
rou_D = zeros(1,length(s_D_E_vec));
for i = 1:length(rou_D)
    rou_D(i) = s_D_I_vec(i)/s_D_E_vec(i);
end

plot(rou_D,ER_vec,'o')
xlim([0 1])

% eigen-analysis

S_mx = S_opt_fit;

[eigen_vec,eigen_val] = eig(S_mx);

% first_eig_val_vec
first_eig_val_vec = [eigen_val_loop1(1,1),eigen_val_loop2(1,1),eigen_val_loop3(1,1),eigen_val_loop4(1,1),eigen_val_loop5(1,1),eigen_val_loop6(1,1),eigen_val_loop21(1,1),eigen_val_loop42(1,1),eigen_val_opt(1,1)];
plot(first_eig_val_vec,s_D_E_vec,'o')

% -second_eig_vec(1)/second_eig_vec(2)

r1 = -eigen_vec_loop1(1,2)/eigen_vec_loop1(2,2);
r2 = -eigen_vec_loop2(1,2)/eigen_vec_loop2(2,2);
r3 = -eigen_vec_loop3(1,2)/eigen_vec_loop3(2,2);
r4 = -eigen_vec_loop4(1,2)/eigen_vec_loop4(2,2);
r5 = -eigen_vec_loop5(1,2)/eigen_vec_loop5(2,2);
r6 = -eigen_vec_loop6(1,2)/eigen_vec_loop6(2,2);
r21 = -eigen_vec_loop21(1,2)/eigen_vec_loop21(2,2);
r42 = -eigen_vec_loop42(1,2)/eigen_vec_loop42(2,2);
ropt = -eigen_vec_opt(1,2)/eigen_vec_opt(2,2);

ratio_second_eig_vec = [r1,r2,r3,r4,r5,r6,r21,r42,ropt];

plot(ratio_second_eig_vec,rou_D,'o')
xlim([0 1])
ylim([0 1])

% eigenval_2_over_S_D_E:

e1 = eigen_val_loop1(2,2)/s_D_E_vec(1);
e2 = eigen_val_loop2(2,2)/s_D_E_vec(2);
e3 = eigen_val_loop3(2,2)/s_D_E_vec(3);
e4 = eigen_val_loop4(2,2)/s_D_E_vec(4);
e5 = eigen_val_loop5(2,2)/s_D_E_vec(5);
e6 = eigen_val_loop6(2,2)/s_D_E_vec(6);
e21 = eigen_val_loop21(2,2)/s_D_E_vec(7);
e42 = eigen_val_loop42(2,2)/s_D_E_vec(8);
eo = eigen_val_opt(2,2)/s_D_E_vec(9);

eigenval_2_over_S_D_E = [e1,e2,e3,e4,e5,e6,e21,e42,eo];



