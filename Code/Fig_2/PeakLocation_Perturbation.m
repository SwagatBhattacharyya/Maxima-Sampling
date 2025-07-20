%% Initialize
clear, clc;

%% Define Constants
w0 = 1; % Keep this 1. This sets the normalized time scale.
a = 1;
P_k = logspace(-4,log(2*pi),1000);
P_tilde_k = P_k/(2*pi);

%% Prepare Signals
delta_tk = zeros(size(P_k));
parfor i = 1:length(P_k)
    Fineness = log(1+abs(P_k(i)))/1e3; % Timestep Precision
    [~,ModulatedSignal,~,t] = PrepareSignals(w0,Fineness,P_k(i),a);
    [~,locs] = findpeaks(ModulatedSignal);
    delta_tk(i) = t(locs);
end

%% Plot
figure;
loglog(P_tilde_k,abs(delta_tk),'r--');
hold on;
loglog(P_tilde_k,abs(P_tilde_k),'k');
loglog(P_tilde_k,abs(P_tilde_k-P_tilde_k.^3),'k-.');
loglog(P_tilde_k,abs(P_tilde_k-P_tilde_k.^3+2*P_tilde_k.^5),'k--');
legend('Numerical Simulation','Analytic: 1^{st} Order','Analytic: 3^{rd} Order','Analytic: 5^{th} Order','location','southeast');
ylabel('|\Deltat_k\omega_0|')
xlabel('$\tilde{P}_k$', 'Interpreter', 'latex');
axis tight
xlim([1e-4 1]);
ylim([1e-4 1]);
xticks(10.^(-4:1:0));
yticks(10.^(-4:1:0));

%% Function to Prepare Signal
function [ModulatingSignal,ModulatedSignal,TheCarrier,t] = PrepareSignals(w0,Fineness,P_k,a) % P_k = percent change in amplitude between carriers, a = starting amplitude
    T = 2*pi/w0;
    t = (-T/2):Fineness:(T/2);
    ModulatingSignal = a+(a*P_k/T)*t;
    TheCarrier = cos(w0.*t);
    ModulatedSignal  = ModulatingSignal.*TheCarrier;
end