%----------------------------------------------- 
% 1-D Zero Range Particle dynamics
% in a Regularized Wiener Random Media.
%
% Euler-Maruyama simulation.
%-----------------------------------------------
clf;

% Perfect seed (periodic traps)
rng(1)

%% Medium
Xstart = 0; Xend = 1; N = 100000;
Xrange = Xstart - Xend;
epsilon = Xrange/N; % regularizing width

% LPSX medium on the Torus = [0, 1].
W = WPmedium(N);

%W.eval(1);
%disp(W.eval(1)*2);
%W.eval(2);
%W.eval(0.2);
%W.eval(-0.2);
%W.eval(-0.3);
%W.eval(-0.17);
%W.eval(-0.173);
%W.eval(-0.175);
%disp(W.eval(0.8) - W.eval(1));
%W.eval(0.8);

plot(W.get_sampling_times(), W.get_path(), 'r-','LineWidth', 2)
xlabel('x_i', 'FontWeight', 'normal')
ylabel('W_i', 'FontWeight', 'normal')
grid on

% Showing extension
figure;
plot([-Xrange:epsilon:Xrange], W.extend(), 'r-', 'LineWidth', 2)
xlabel('x_i', 'FontWeight', 'normal')
ylabel('W_i', 'FontWeight', 'normal')
grid on

%% Speed change
%psi = @(t,x) 3 + sin(x) + cos(t);

%% Euler-Maruyama
Xzero = 0;
L = 10000; T = 100; dt = T/L;

Xem = zeros(1,L+1);
Xem(1) = Xzero;
for j = 1:L
    dBt = sqrt(dt) * randn; % dBt

    % 1-D Brox with LPSX medium
    Xem(j+1) = Xem(j) + dBt - W.d_dx(Xem(j), epsilon) * dt;
end

figure;
plot([0:dt:T], Xem, 'k-', 'LineWidth', 1)
xlabel('t', FontSize=12)
ylabel('Xem', 'FontWeight', 'normal')
grid on