%% Initialize MATLAB
clear, clc;

%% Define Constants
w_0 = 1;
Q = 10;
a = w_0/(2*Q);
t = 0:1e-6:(6*pi+0.001); % Three Periods

%% Compute Curves
M = exp(-a*t); % Compute Modulating Function
V_in = cos(w_0*t).*M;
Base = exp(-2*pi*a/w_0);
Levels = Base.^[0 1 2 3];
StepCurve = zeros(size(t));
for i = 1:length(Levels)
StepCurve((t >= 2*pi*(i-1)) & (t < 2*pi*i)) = Levels(i);
end

%% Plot
figure, hold on;
plot(t,V_in,'k');
plot(t,M,'k--')
plot(t,StepCurve,'b')
ylim([-0.9 1.1]);
xlim([0 6*pi+0.2]);
axis tight;