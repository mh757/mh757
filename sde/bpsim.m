% Brownian Path Simulation
% ----------------------------------------------
% Programming Exercise in Chapter 03 and 04 in
% An Introduction to the Numerical Simulation of
% Stochastic Differential Equations
% by Higham and Kloeden.
% ----------------------------------------------
% Author: Marcel Hudiani
% 10/27/2024
% ----------------------------------------------

%% Chapter 3
clf
rng(100)
L = 500; T = 1; dt = T/L;
W = WP(0, L, T); % BM starting at the origin.
W.print();

plot(W.get_sampling_times(), W.get_path(), 'r-','LineWidth', 2)
xlabel('t_i', 'FontWeight', 'normal')
ylabel('W_i', 'FontWeight', 'normal')
grid on
hold on

% Now, increase the sample frequency by 2.
W = W.fill();
W.print(); % Print out the new object

% Draw the path in the same figure.
plot(W.get_sampling_times(), W.get_path(), 'r--o','LineWidth', 2, ...
    'Color', [0 1 0], 'MarkerSize', 4, 'MarkerFaceColor', 'auto')
hold off


%% Chapter 4
% Compute the Ito integral of W dW
integrand = W.get_path();
[ito, itoerr] = W.integrate(integrand(1 : end - 1), false);

disp(['The Ito integral is ', num2str(ito), ...
    ' and the error is ', num2str(itoerr)])