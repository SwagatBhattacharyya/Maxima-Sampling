%% Initialize
clear, clc;

%% Define Constants
tspan = [0,1.5];
tstart = 0.1;
tend = 0.6;
Params.tau_att = 0.001;
TrackingLevel = 100; % High Ripple
%TrackingLevel = 500; % Low Ripple
Params.tau_dec = Params.tau_att*TrackingLevel;

%% Generate Signals
Params.t = tspan(1):1e-6:tspan(2);
Params.Vin = cos((2*pi/0.05)*(Params.t-tstart));
Params.Vin(Params.t < tstart) = 0;
Params.Vin(Params.t > tend) = 0;

%% Simulate
y0 = 0;
[t,y] = ode15s(@(t,y)PKDE(t,y,Params),Params.t,y0);

%% Plot
figure, hold on;
plot(Params.t,Params.Vin,'k');
plot(t,y,'r');
axis tight;

%% ODE Function
function dydt = PKDE(t,y,Params)
% Interpolate the input
Vin = interp1(Params.t,Params.Vin,t,"linear");

% Compute update
if(Vin > y)
    tau = Params.tau_att;
else
    tau = Params.tau_dec;
end
dydt = (Vin-y)/tau;
end