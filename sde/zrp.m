%-----------------------------------------------
% 1-D ZRP in Regularized, Wiener Random Media 
% 1-D Particle dynamics
%-----------------------------------------------

% 1-D Random medium (spatial)
W = bm(0, 1, 'StartTime', -1000);

% Set seed and Random Number Generator
rng(10203,'twister')

% Simulate W to obtain realization so that we may treat them as functions.
dx       = 1 / 10;
nPeriods = 20000;
X        = nPeriods * dx;
nPaths   = 50;    % # of simulated paths (50 data points for ML)
sampleTimes = cumsum([W.StartTime; dx(ones(nPeriods,1))]);
Wx = W.simulate(nPeriods, 'DeltaTime', dx, 'nTrials', nPaths, 'Z', z);

%TODO: What's z? z = Example_StratifiedRNG(nPaths, sampleTimes);


%--------------------------------------
% Particle 
%--------------------------------------

% Regularization of the medium
eReg = 0.01;
W_reg = (W(x+eReg) - W(x-eReg))/eReg

% Speed change
psi = @(t,x) 3 + sin(x) + cos(t);

% Particle dynamics is an SDE:
% dxt = dBt - 1/2 W_reg(xt) * psi dt
F = @(t,x) -0.5 * psi;
G = @(t,x) 1;

xt_sde = sde(F, G);