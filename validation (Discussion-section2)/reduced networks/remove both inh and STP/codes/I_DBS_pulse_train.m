%% ===============================
%  DBS Pulse Train Generator
%  ===============================

% clear; clc;

%% ===== User-defined parameters =====
f_DBS      = 130;      % DBS frequency (Hz)
T_total    = 10;       % total stimulation time (s)
dt         = 0.1e-3;   % sampling resolution (s) -> 0.1 ms
pulse_w    = 0.3e-3;   % pulse width (s) -> 0.3 ms
A          = 1;        % pulse amplitude (scalar)

%% ===== Derived parameters =====
Fs         = 1/dt;                 % sampling frequency
N_total    = round(T_total/dt);    % total number of samples
t          = (0:N_total-1) * dt;   % time vector

T_period   = 1/f_DBS;              % period of DBS
N_period   = round(T_period/dt);   % samples per period

N_pulse    = round(pulse_w/dt);    % samples per pulse
% (0.3 ms / 0.1 ms = 3 samples; 
% if you want exactly 4 samples, set N_pulse = 4 manually)

%% ===== Generate pulse train =====
dbs = zeros(1, N_total);

for k = 1:N_period:N_total
    idx_end = min(k + N_pulse - 1, N_total);
    dbs(k:idx_end) = A;
end

%% ===== Plot first 100 ms for visualization =====
figure;
plot(t*1000, dbs, 'LineWidth', 1.5);
xlim([0 100]);  % show first 100 ms
xlabel('Time (ms)');
ylabel('Amplitude');
title(['DBS Pulse Train (' num2str(f_DBS) ' Hz)']);
set(gca,'FontSize',12,'LineWidth',1.2);