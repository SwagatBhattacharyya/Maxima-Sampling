%% Creates a random bandlimited synthetic signal and performs Max-SED analysis on it
%% Init
clear,  clc;

%% Define Constants
ModIdx = 0.95;      % Modulation Index
w0 = 1;             % Keep this 1. This sets the normalized time scale.
Fs = 1e3*w0/(2*pi); % Omega_sample = 1000x carrier freq for convenience 
N = 2^22;           % N point DFT
l = 100;            % l < N/2 sets omega_l (bandlimit)
H = 0.01;           % Random number generation range 

%% Compute signals
[omega_l, UkEst, M, t] = ConstuctRandom(l, N, Fs, ModIdx, w0,H);
ModulatedSignal = M.*cos(w0.*t);
Demodulated = DemodulateArbitrary(ModulatedSignal);
NormalizedPower = ErrorPower(M,Demodulated,ModulatedSignal);
max(NormalizedPower)

%% Plot
figure, hold on;
plot(t,Demodulated);
plot(t,M,'r--');
plot(t,ModulatedSignal,'k');
legend('V_{out}','M(t)','V_{in}=M(t)cos(\omega_0t)');

% Formatting
axis tight;
xlabel('t (s)')
ylabel('Volts')

%% Constuct Signal Fcn
function [omega_l, U_k_est, m, t] = ConstuctRandom(l, N, Fs, ModIdx, w0,H)
    %% Allocate full spectrum (row-vector)
    M_freq          = zeros(1,N);            % k = 0 … N–1

    % k = 1 … l  (positive-frequency bins <= ω_l)
    posBins         = randn(1,l-1) + 1j*randn(1,l-1);
    M_freq(2:l)     = N * posBins;           % scale for ifft convention

    % Negative-frequency replicas: k = N-l+2 … N
    M_freq(N-l+2:N) = conj(M_freq(l:-1:2));

    %% Time-domain signal
    m   = ifft(M_freq,'symmetric');          % length N
    t   = (0:N-1)/Fs;                        % seconds

    %% Stretch and scale
    m   = m - mean([min(m),max(m)]);
    m = 1+2*ModIdx*m/range(m);

    %% Spectral figures
    omega_l = 2*pi*l*Fs/N;               % rad/s

    %% error bound
    U_k_est = (4*(pi^2)/3)*(omega_l/w0)^2 *(log(1-ModIdx)).^2;
end

%% Demodulate
function DemodulatedSignal = DemodulateArbitrary(ModulatedSignal)
    % Find the peaks
    [pks,idx] = findpeaks(ModulatedSignal);

    % Simulate sample and hold
    DemodulatedSignal = interp1(idx,pks,1:length(ModulatedSignal),"previous");
end

%% Define Piecewise Error Computation Function
function NormalizedPower = ErrorPower(Original,Demodulated,ModulatedSignal)
    [~,pkidx] = findpeaks(ModulatedSignal);
    NormalizedPower = zeros([1,length(pkidx)-1]);
    for i = 1:(length(pkidx)-1)
        MaskIdx = pkidx(i):pkidx(i+1);
        Power = mean((Original(MaskIdx)-Demodulated(MaskIdx)).^2);
        Den = mean(Demodulated(MaskIdx).^2);
        NormalizedPower(i) = Power/Den;
    end
end