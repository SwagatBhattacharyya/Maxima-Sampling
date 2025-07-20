%% Creates a random bandlimited synthetic signal and performs Max-SED analysis on it
%% Init
clear,  clc;

%% Define Constants
w0 = 1;                     % Keep this 1. This sets the normalized time scale.
Fs = 1.668860536075273e+02; % Omega_sample = roughly 1000x carrier freq for convenience 
N = 2^20;                   % N point DFT

%% Define Sweep Range
All_l = [10,100,1000];
All_ModIdx = logspace(-3,log10(0.9),25);
Trials = 24;

%% Compute signals and test
nL   = numel(All_l);
nM   = numel(All_ModIdx);
Bound     = zeros(nL,nM,Trials);
Numerical = zeros(nL,nM,Trials);

parfor k = 1:Trials
    for i = 1:nL
        for j = 1:nM
            [~,UkEst,M,t] = ConstuctRandom(All_l(i),N,Fs,All_ModIdx(j),w0);
            X = M.*cos(w0*t);
            Y = DemodulateArbitrary(X);
            Bound(i,j,k)     = UkEst;
            Numerical(i,j,k) = max(ErrorPower(M,Y,X));
        end
    end
end

%% Compute W
for i = 1:nL
    [myB,UkEst,M,t] = ConstuctRandom(All_l(i),N,Fs,All_ModIdx(1),w0);
    w(i) = myB/w0;
end

%% Plot
colors  = {'r','b','k'};         % colour per l (i-index)
figure;
ax = gca;
set(ax,'XScale','log','YScale','log');   % both axes logarithmic
hold(ax,'on');

% Bounds
for i = fliplr(1:nL)
    yLine = squeeze(Bound(i,:,1));                 % size 1 × nM
    plot(All_ModIdx, yLine, [colors{i} '--'], ...
         'LineWidth', 1);
end
for i = 1:nL
    for j = 1:nM
        yVals = squeeze(Numerical(i,j,:));              % Trials × 1
        xVals = All_ModIdx(j) .* ones(size(yVals));     % same x for all k

        % draw the box
        h = boxchart(xVals, yVals, ...
                     'BoxFaceColor', colors{i}, ...
                     'WhiskerLineColor', colors{i}, ...
                     'MarkerStyle', 'none');
        % thinner boxes improve visibility on log axes btw
        h.BoxWidth = 0.25*All_ModIdx(j);   % width ≈ 4 % of x-location
    end
end

% Formatting
axis tight
xlabel('m')
ylabel('max(U_k)')
lgd = arrayfun(@(idx) ...
       sprintf('$w = %d$', w(idx)), ...
       fliplr(1:nL), 'UniformOutput', false);
legend('w=1.0','w=0.1','w=0.01', 'Interpreter','latex', 'Location','best');
yticks(10.^(-10:2:2))
xlim([1e-3,1])
ylim([10^-10,10^2])

%% Constuct Signal Fcn
function [omega_l, U_k_est, m, t] = ConstuctRandom(l, N, Fs, ModIdx, w0)
    %% Allocate full spectrum
    M_freq          = zeros(1,N);            % k = 0 … N–1

    % k = 1 … l  (positive-frequency bins <= omega_l)
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

    %% Error bound
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