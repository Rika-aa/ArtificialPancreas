clc
clear
close all

%% MODEL

%Type 1 Diabtes Model Parameters
p1 = 0.0287;
p2 = 0.0283;
p3 = 5.035*10^(-5);
p4 = 5/54;

Gb = 119.1858;
Ib = 15.3872;


%% EQUILIBRIUM

%Initial Conditions
x0 = [Gb, 0 ,Ib]';
u0 = [0; 0];
y0 = x0;
idu = [1; 2];

%Equilibrium point
[xeq, ueq] = trim('Pianta', x0, u0, y0, [] ,idu,[])
yeq = xeq(1,:);

%% LINEARIZATION

%linearization of the model
[A,B,C,D] = linmod('Pianta', xeq,ueq)

%C = C(1,:)
%D = D(1,:)

%State Space Representation
sys = ss(A,B,C,D);


%% SYSTEM ANALYSIS
%Eigenvalues of matrix A
eigenvalues = eig(A);
[V,E,W] = eig(A);

% %% Plot risposta libera sys lin
%
% out = sim('Modello_Lineare');
% t = out.xout.Time;                 
% 
% x_data = squeeze(out.xout.Data)'; 
% x1 = x_data(:, 1);          
% x2 = x_data(:, 2);         
% x3 = x_data(:, 3);        
% figure('Name', 'Free Response of the System - Case 1', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
%
% % x1(t) 
% subplot(3, 1, 1);
% plot(t, x1, 'b', 'LineWidth', 1.5);
% title('Free Response: Glucose Variation', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('$\tilde{x}_1(t)$ [mg/dL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %x2(t)
% subplot(3, 1, 2);
% plot(t, x2, 'r', 'LineWidth', 1.5);
% title('Free Response: Active Insulin Variation', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('$\tilde{x}_2(t)$ [1/min]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %x3(t)
% subplot(3, 1, 3);
% plot(t, x3, 'g', 'LineWidth', 1.5);
% title('Free Response: Plasma Insulin Variation', 'Interpreter', 'latex', 'FontSize', 12);
% xlabel('Time [min]', 'Interpreter', 'latex', 'FontSize', 11);
% ylabel('$\tilde{x}_3(t)$ [$\mu$U/mL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;

%% FORCED RESPONSE (1)
Bu = [0; 0; 1]; 
Bd = [1; 0; 0]; 
C = [1, 0, 0]; 
Du = 0;
Dd = 0;

%State-Space
sys_u = ss(A, Bu, C, Du); 
sys_d = ss(A, Bd, C, Dd);

%Funzioni di Trasferimento
Wu = tf(sys_u);
Wd = tf(sys_d);

% %Plot Forced Response
% T_sim = 300; 
% 
% [y_u, t_u] = step(sys_u, T_sim); 
% [y_d, t_d] = step(sys_d, T_sim);
% 
% figure('Name', 'Forced Response Analysis', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
% 
% %u(t)
% subplot(2, 1, 1);
% plot(t_u, y_u, 'b', 'LineWidth', 1.5);
% title('Forced Response to Insulin Input $u(t)$ [$d(t) = 0$]', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('$\tilde{x}_1(t)$ [mg/dL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %d(t)
% subplot(2, 1, 2);
% plot(t_d, y_d, 'r', 'LineWidth', 1.5);
% title('Forced Response to Meal Disturbance $d(t)$ [$u(t) = 0$]', 'Interpreter', 'latex', 'FontSize', 12);
% xlabel('T [min]', 'Interpreter', 'latex', 'FontSize', 11);
% ylabel('$\tilde{x}_1(t)$ [mg/dL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;

%Poli
Pu = pole(Wu);
Pd = pole(Wd);


%% %FORCED RESPONSE (2)
% iu_1 = 1;
% iu_2 = 2;
% 
% [num1 den1] = ss2tf(A,B,C,D(1,:),iu_1);
% [num2 den2] = ss2tf(A,B,C,D(1,:),iu_2);
% 
% W_yu_1 = tf(num1,den1)
% W_yu_2 = tf(num2,den2)


%% STRUCURAL ANALYSIS
R = ctrb(A,Bu);
rank_R = rank(R);

%% FULL STATE FEEDBACK
%poles_SSF = [-0.2, -0.3, -0.4];
poles_SSF = [-0.05, -0.06, -0.1];
K = place(A, Bu, poles_SSF);
eig(A-Bu*K);

% Open Loop (Paziente diabetico)
% modelSim_OFF = sim('Sys_OL');
% time_off = modelSim_OFF.tout;
% x1_off = modelSim_OFF.x_out(:,1); %G(t)
% x2_off = modelSim_OFF.x_out(:,2); %X(t)
% x3_off = modelSim_OFF.x_out(:,3); %I(t)
% 
% % Closed Loop (Pancreas Artificiale)
% modelSim_ON = sim('Sys_CL');
% time_on = modelSim_ON.tout;
% x1_on = modelSim_ON.x_out(:,1); %G(t)
% x2_on = modelSim_ON.x_out(:,2); %X(t)
% x3_on = modelSim_ON.x_out(:,3); %I(t)
% 
% %Confronto G(t)
% figure;
% plot(time_off, x1_off, 'r--', 'LineWidth', 1.5);
% hold on;
% plot(time_on, x1_on, 'b', 'LineWidth', 1.5);
% hold off;
% title('Comparison of Glucose concentration $G(t)$','Interpreter','latex')
% legend('Open Loop', 'Closed Loop ($u = K\tilde{x}$)','Interpreter', 'latex', 'Location', 'best');
% xlabel('T (min)');
% ylabel('$G(t)$', 'Interpreter', 'latex');
% grid on;
% 
% %Confronto X(t)
% figure;
% plot(time_off, x2_off, 'r--', 'LineWidth', 1.5);
% hold on;
% plot(time_on, x2_on, 'b', 'LineWidth', 1.5);
% hold off;
% title('Comparison of Insulin Action $X(t)$','Interpreter','latex')
% legend('Open Loop', 'Closed Loop ($u = K\tilde{x}$)','Interpreter', 'latex', 'Location', 'best');
% xlabel('T (min)');
% ylabel('$X(t)$', 'Interpreter', 'latex');
% grid on;
% 
% %Confronto I(t)
% figure;
% plot(time_off, x3_off, 'r--', 'LineWidth', 1.5);
% hold on;
% plot(time_on, x3_on, 'b', 'LineWidth', 1.5);
% hold off;
% title('Comparison of Plasma Insulin $I(t)$','Interpreter','latex')
% legend('Open Loop', 'Closed Loop ($u = K\tilde{x}$)','Interpreter', 'latex', 'Location', 'best');
% xlabel('T (min)');
% ylabel('$I(t)$', 'Interpreter', 'latex');
% grid on;


%% LUENBERGER
O = obsv(A,C);
RankO = rank(O);

%poles_obs = 2* [-0.2, -0.3, -0.4];
poles_obs = [-0.15, -0.20, -0.30];

Lt = place(A',C',poles_obs);
L=Lt';
AL = A-L*C;
eig(AL);
BL = [Bu L];
CL = eye(size(A));
DL = zeros(size(A,1),2);


% % Plot
% out = sim('Obs_L');
% time = out.tout;
% 
% x_real_tilde_1 = out.x_real_tilde(:,1);
% x_real_tilde_2 = out.x_real_tilde(:,2);
% x_real_tilde_3 = out.x_real_tilde(:,3);
% 
% x_est_tilde_1 = out.x_est_tilde(:,1);
% x_est_tilde_2 = out.x_est_tilde(:,2);
% x_est_tilde_3 = out.x_est_tilde(:,3);
% 
% 
% error_1 = x_real_tilde_1 - x_est_tilde_1;
% error_2 = x_real_tilde_2 - x_est_tilde_2;
% error_3 = x_real_tilde_3 - x_est_tilde_3;
% 
% %Variazioni Glicemia
% figure;
% plot(time, x_real_tilde_1, 'k', 'LineWidth', 1.5);
% hold on;
% plot(time, x_est_tilde_1, 'b--', 'LineWidth', 1.5);
% hold off;
% title('Variation Comparison: Real vs Estimated Glucose $\tilde{x}_1(t)$', 'Interpreter', 'latex');
% legend('Real variation $\tilde{x}_1(t)$', 'Estimated variation $\hat{\tilde{x}}_1(t)$', ...
%        'Interpreter', 'latex', 'Location', 'best');
% xlabel('time [min]');
% ylabel('$\tilde{x}_1(t)$ e $\hat{x}_1(t)$ [mg/dL]', 'Interpreter', 'latex');
% grid on;
% 
% %Variazioni Insulina Attiva
% figure;
% plot(time, x_real_tilde_2, 'k', 'LineWidth', 1.5);
% hold on;
% plot(time, x_est_tilde_2, 'b--', 'LineWidth', 1.5);
% hold off;
% title('Variation Comparison: Real vs Estimated Active Insulin $\tilde{x}_2(t)$', 'Interpreter', 'latex');
% legend('Real variation $\tilde{x}_2(t)$', 'Estimated variation $\hat{x}_2(t)$', ...
%        'Interpreter', 'latex', 'Location', 'best');
% xlabel('time [min]');
% ylabel('$\tilde{x}_2(t)$ e $\hat{x}_2(t)$ [1/min]', 'Interpreter', 'latex');
% grid on;
% 
% %Variazioni: Plasmatica
% figure;
% plot(time, x_real_tilde_3, 'k', 'LineWidth', 1.5);
% hold on;
% plot(time, x_est_tilde_3, 'b--', 'LineWidth', 1.5);
% hold off;
% title('Variation Comparison: Real vs Estimated Plasma Insulin $\tilde{x}_3(t)$', 'Interpreter', 'latex');
% legend('Real variation $\tilde{x}_3(t)$', 'Estimated variation $\hat{x}_3(t)$', ...
%        'Interpreter', 'latex', 'Location', 'best');
% xlabel('time [min]');
% ylabel('$\tilde{x}_3(t)$ e $\hat{x}_3(t)$ [mU/L]', 'Interpreter', 'latex');
% grid on;
% 
% %Errore di Stima: e_1(t)
% figure;
% plot(time, error_1, 'r', 'LineWidth', 1.5);
% title('Estimation Error for Glucose: $e_1(t) = \tilde{x}_1(t) - \hat{x}_1(t)$', 'Interpreter', 'latex');
% xlabel('time [min]');
% ylabel('$e_1(t)$', 'Interpreter', 'latex');
% grid on;
% 
% %Errore di Stima: e_2(t)
% figure;
% plot(time, error_2, 'r', 'LineWidth', 1.5);
% title('Estimation Error for Insulin Action: $e_2(t) = \tilde{x}_2(t) - \hat{x}_2(t)$', 'Interpreter', 'latex');
% xlabel('time [min]');
% ylabel('$e_2(t)$', 'Interpreter', 'latex');
% grid on;
% 
% %Errore di Stima: e_3(t)
% figure;
% plot(time, error_3, 'r', 'LineWidth', 1.5);
% title('Estimation Error for Plasma Insulin: $e_3(t) = \tilde{x}_3(t) - \hat{x}_3(t)$', 'Interpreter', 'latex');
% xlabel('time [min]');
% ylabel('$e_3(t)$', 'Interpreter', 'latex');
% grid on;


%% Dynamic Compensator
% out = sim('Dynamic_Compensator');
% x1d = out.xout.Data(:,1); 
% x2d = out.xout.Data(:,2); 
% x3d = out.xout.Data(:,3);
% time_d = out.xout.Time;
% 
% figure('Name', 'Dynamic Compensator (non-linear)', 'NumberTitle', 'off', 'Position', [100, 100, 800, 800]);
% 
% %x1
% subplot(3, 1, 1);
% plot(time_d, x1d, 'b', 'LineWidth', 1.5);
% title('Dynamic Compensator: Glucose Variation $x_1(t)$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('[mg/dL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %x2
% subplot(3, 1, 2);
% plot(time_d, x2d, 'r', 'LineWidth', 1.5);
% title('Dynamic Compensator: Active Insulin Variation $x_2(t)$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('[1/min]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %x3
% subplot(3, 1, 3);
% plot(time_d, x3d, 'g', 'LineWidth', 1.5);
% title('Dynamic Compensator: Plasma Insulin Variation $x_3(t)$', 'Interpreter', 'latex', 'FontSize', 12);
% xlabel('T [min]', 'Interpreter', 'latex', 'FontSize', 11);
% ylabel('[$\mu$U/mL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %% Dynamic Compensator - lineare
% out = sim('Dynamic_Compensator_Linear');
% x1d = out.xout.Data(:,1); 
% x2d = out.xout.Data(:,2); 
% x3d = out.xout.Data(:,3);
% time_d = out.xout.Time;
% 
% figure('Name', 'Dynamic Compensator (linear)', 'NumberTitle', 'off', 'Position', [100, 100, 800, 800]);
% 
% %x1
% subplot(3, 1, 1);
% plot(time_d, x1d, 'b', 'LineWidth', 1.5);
% title('Dynamic Compensator: Glucose Variation $x_1(t)$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('[mg/dL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% % x2
% subplot(3, 1, 2);
% plot(time_d, x2d, 'r', 'LineWidth', 1.5);
% title('Dynamic Compensator: Active Insulin Variation $x_2(t)$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('[1/min]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;
% 
% %x3
% subplot(3, 1, 3);
% plot(time_d, x3d, 'g', 'LineWidth', 1.5);
% title('Dynamic Compensator: Plasma Insulin Variation $x_3(t)$', 'Interpreter', 'latex', 'FontSize', 12);
% xlabel('T [min]', 'Interpreter', 'latex', 'FontSize', 11);
% ylabel('[$\mu$U/mL]', 'Interpreter', 'latex', 'FontSize', 11);
% grid on;