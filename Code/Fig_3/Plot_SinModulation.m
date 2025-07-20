%% Initialize
clear, clc;

%% Plot
figure;

% 0.1
load('wm0_1.mat');
Expected_SecondOrder = (4*pi^2/3).*P_tilde_k.^2;
loglog(m,NormalizedErrorPowers,'m--')
hold on;
loglog(m,Expected_SecondOrder,'m')

% 0.01
load('wm0_01.mat');
Expected_SecondOrder = (4*pi^2/3).*P_tilde_k.^2;
loglog(m,NormalizedErrorPowers,'b--')
loglog(m,Expected_SecondOrder,'b')

% 0.001
load('wm0_001.mat');
Expected_SecondOrder = (4*pi^2/3).*P_tilde_k.^2;
loglog(m,NormalizedErrorPowers,'k--')
loglog(m,Expected_SecondOrder,'k')

% Formatting
legend('Sim: W=0.1','Analytic: W=0.1','Sim: w=0.01','Analytic: w=0.01','Sim: w=0.001','Analytic: w=0.001','location','northwest')
ylabel('max$(U_k)$', 'Interpreter', 'latex')
xlabel('$m$', 'Interpreter', 'latex');
axis tight
xlim([1e-4 1]);
ylim([1e-13 1]);
xticks(10.^(-4:1:0));
yticks(10.^(-14:2:1));