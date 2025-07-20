%% This uses numerical integration to compute how the demodulation performance of a maxima sampling envelope detector varies as a modulation index
%% Initialize
clear, clc;

%% Define Simulation Hyperparameters
wm = 0.001;
FileName = 'wm0_001.mat';

%% Define Constants
w0 = 1; % Keep this 1. This sets the normalized time scale.
Fineness = 1e-4; % Integration precision
m = [logspace(-6,-1,250),logspace(-1,0,250)];
P_tilde_k = (wm./w0).*m./sqrt(1-m.^2);

%% Run parametric sweep
a = 1;
NormalizedErrorPowers = zeros(size(P_tilde_k));
parfor i = 1:length(P_tilde_k)
    [ModulatingSignal,ModulatedSignal,t] = PrepareSignals(w0,Fineness,m(i),wm);
    Demodulated = DemodulateArbitrary(ModulatedSignal);
    NormalizedPower = ErrorPower(ModulatingSignal,Demodulated,ModulatedSignal);
    NormalizedErrorPowers(i) = max(NormalizedPower);
end
save(FileName);

%% Plot
Expected_SecondOrder = (4*pi^2/3).*P_tilde_k.^2;
figure;
loglog(m,NormalizedErrorPowers,'r--')
hold on;
loglog(m,Expected_SecondOrder,'k')
loglog([1e-3,1e-2],[1e-4,1e-2],'k') % 20 dec/dec
legend('Numerical Simulation','Analytic: 2^{nd} Order','location','southeast')
ylabel('max$(U_k)$', 'Interpreter', 'latex')
xlabel('$m$', 'Interpreter', 'latex');
axis tight

%% Demodulate
function DemodulatedSignal = DemodulateArbitrary(ModulatedSignal)
    % Find the peaks
    [pks,idx] = findpeaks(ModulatedSignal);

    % simulate sample and hold
    DemodulatedSignal = interp1(idx,pks,1:length(ModulatedSignal),"previous");
end

%% Define Piecewise Error Computation Function
function NormalizedPower = ErrorPower(Original,Demodulated,ModulatedSignal)
    % Compute
    [~,pkidx] = findpeaks(ModulatedSignal);
    NormalizedPower = zeros([1,length(pkidx)-1]);
    for i = 1:(length(pkidx)-1)
        MaskIdx = pkidx(i):pkidx(i+1);
        Power = mean((Original(MaskIdx)-Demodulated(MaskIdx)).^2);
        Den = mean(Demodulated(MaskIdx).^2);
        NormalizedPower(i) = Power/Den;
    end
end

%% Prepare Signals
function [ModulatingSignal,ModulatedSignal,t] = PrepareSignals(w0,Fineness,m,wm)
    T = 2*pi/wm;
    t = 0:Fineness:T;
    ModulatingSignal = 1 + m*cos(wm.*t);
    ModulatedSignal  = ModulatingSignal.*cos(w0.*t);
end