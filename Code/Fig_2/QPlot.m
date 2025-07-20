%% Initialize
clear, clc;

%% Plot Solution of Q vs U for full analytic expression
Q = logspace(-2, 4, 1000);
U_Analytical_Time = (4*exp(-pi./Q)-exp(-2*pi./Q)+2*pi./Q - 3)./(1-exp(-2*pi./Q));

%% Compute from numerical integration
SimulationTimeConstants = 100; 
Q_Numerical = logspace(-2,4, 100);
U_Numerical = zeros(size(Q_Numerical));
dt = 1e-5;
parfor i = 1:length(Q_Numerical)
    % Set Up
    myQ = Q_Numerical(i);
    tspan = [0,SimulationTimeConstants*2*myQ];
    a = 1/(2*myQ);

    % Find peaks
    t = (0:dt:SimulationTimeConstants)*2*myQ;
    TrueEnvelope = exp(-t/(2*myQ));
    Vin = cos(t).*TrueEnvelope;
    [peaks,idx] = findpeaks(Vin);
    peaktimestamps = [0,t(idx)];
    peaks = [Vin(1),peaks]; % Approximately.

    % Integrate
    [~,E] = ode89(@(t,y)ODEFun(t,a,peaks),tspan,0)
    U_Numerical(i) = E(end)*(2*a);
end

%% Compute bound from linearization
Short_Q = logspace(-1, 4, 1000);
LinearizationBound = (pi^2)./(3.*(Short_Q.^2)) + pi./(8.*(Short_Q.^3));

%% Plot
figure;
loglog(Q,U_Analytical_Time,'b');         
hold on;
loglog(Q_Numerical,U_Numerical,'r--') % From numerical integration
loglog(Short_Q,LinearizationBound,'k') % From time domain bound
ylabel('U');
xlabel('Q');
legend('Analytical Result', 'Numerical Simulation', 'Linearization-Based Bound','location','southwest')
yticks(10.^(-7:2:3))
xticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4])
xlim([0.01,1e4]);
ylim([1e-7, 1e3]);

%% ODE Function
function dEdt = ODEFun(t,a,peaks)
    %% Extract, assuming Omega = 1;
    T = 2*pi;
    k = floor(t/T);

    %% Compute dEdt
    dEdt = (peaks(k+1)-exp(-a*t)).^2;
end