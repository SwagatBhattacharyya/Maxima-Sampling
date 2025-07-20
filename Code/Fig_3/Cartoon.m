%% Initialize
clear, clc;

%% Define Constants
w0 = 1; % Keep this 1. This sets the normalized time scale.
a = 1;
Fineness = 1e-3; % Timestep Precision
P_k = 1;
P_tilde_k = P_k/(2*pi);

%% Prepare Signals
[ModulatingSignal,ModulatedSignal,TheCarrier,t] = PrepareSignals(w0,Fineness,P_k,a);
[peaks,locs] = findpeaks(ModulatedSignal);
tpeaks = t(locs);
prevpeak = peaks(1)-diff(peaks);

%% Plot
figure;
hold on;
plot(t,ModulatedSignal,'k');
plot(t,ModulatingSignal,'r--');
plot(t,TheCarrier,'k--');
xline([0 2*pi],'k--')
xline(tpeaks,'k')
plot([min(t),tpeaks(1),tpeaks(1)+Fineness,tpeaks(2),tpeaks(2)+Fineness,max(t)],[prevpeak,prevpeak,peaks(1),peaks(1),peaks(2),peaks(2)],'b');

%% Function to Prepare Signal
function [ModulatingSignal,ModulatedSignal,TheCarrier,t] = PrepareSignals(w0,Fineness,P_k,a) % P_k = percent change in amplitude between carriers, a = starting amplitude
    T = 2*pi/w0;
    t = (-T/4):Fineness:(5*T/4);
    ModulatingSignal = a+(a*P_k/T)*t;
    TheCarrier = cos(w0.*t);
    ModulatedSignal  = ModulatingSignal.*TheCarrier;
end