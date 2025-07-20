%% This uses numerical integration to compute how the demodulation performance of a maxima sampling envelope detector varies as a function of the rate of increase between adjacent peaks

%% Initialize
clear, clc;

%% Define Constants
w0 = 1; % Keep this 1. This sets the normalized time scale.
Fineness = 1e-3; % Integration precision
P_k = logspace(-4,log(2*pi),100);
P_tilde_k = P_k/(2*pi);

%% Run parametric sweep
a = 1;
NormalizedErrorPowers = zeros(size(P_k));
for i = 1:length(P_k)
    [ModulatingSignal,ModulatedSignal,t] = PrepareSignals(w0,Fineness,P_k(i),a);
    [pks,idx] = findpeaks(ModulatedSignal);
    Demodulated = pks(1)*ones(size(ModulatingSignal)); % Simple zero order hold
    tspan = t(idx);
    [~,NormalizedPower] = ErrorPower(ModulatingSignal,Demodulated,t,tspan,a);
    NormalizedErrorPowers(i) = NormalizedPower;  
end

%% Plot
Expected_FourthOrder = (4*pi^2/3).*P_tilde_k.^2 + pi*(P_tilde_k.^3) + (1/4 - 8*pi^2/3).*(P_tilde_k.^4);
Expected_ThirdOrder = (4*pi^2/3).*P_tilde_k.^2 + pi*(P_tilde_k.^3);
Expected_SecondOrder = (4*pi^2/3).*P_tilde_k.^2;
figure;
loglog(P_tilde_k,NormalizedErrorPowers,'r--')
hold on;
loglog(P_tilde_k,Expected_SecondOrder,'k')
loglog(P_tilde_k,Expected_ThirdOrder,'k-.')
loglog(P_tilde_k,Expected_FourthOrder,'k--')
legend('Numerical Simulation','Analytic: 2^{nd} Order','Analytic: 3^{rd} Order','Analytic: 4^{th} Order','location','southeast')
ylabel('$U_k$', 'Interpreter', 'latex')
xlabel('$\tilde{P}_k$', 'Interpreter', 'latex');
axis tight
xlim([1e-4 1]);
ylim([1e-7 10]);
xticks(10.^(-4:1:0));
yticks(10.^(-7:2:1));

%% Define Error Computation Function
function [Power,NormalizedPower] = ErrorPower(Original,Demodulated,t,tspan,a)
    MaskIdx = (t>tspan(1)) & (t<tspan(2));
    Power = mean((Original(MaskIdx)-Demodulated(MaskIdx)).^2);
    Den = mean(Demodulated(MaskIdx).^2);
    NormalizedPower = Power/Den;
end

%% Prepare Signals
function [ModulatingSignal,ModulatedSignal,t] = PrepareSignals(w0,Fineness,P_k,a) % P_k = percent change in amplitude between carriers, a = starting amplitude
    T = 2*pi/w0;
    t = (-T/2):Fineness:(3*T/2);
    ModulatingSignal = a+(a*P_k/T)*t;
    ModulatedSignal  = ModulatingSignal.*cos(w0.*t);
end