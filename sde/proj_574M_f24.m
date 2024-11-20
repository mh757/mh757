%----------------------------------------------- 
% 1-D Zero Range Particle dynamics
% in a Regularized Wiener Random Media.
%
% Euler-Maruyama simulation.
%-----------------------------------------------
clf;

% DEBUG
ISPLOT = false;
ISWRITE = true;

% Perfect seed (periodic traps)
seednum = 1;
rng(seednum)

%% Medium
Xstart = 0; Xend = 1; N = 100000;
Xrange = Xstart - Xend;
epsilon = Xrange/N; % regularizing width

% LPSX medium on the Torus = [0, 1].
W = WPmedium(N);

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
psi = @(t,x,m) m + sin(x) + cos(t);

%% Euler-Maruyama
Xzero = 0;
L = 200; T = 10; dt = T/L;
t = 0:dt:T;

% BM
B = WP(0,200,10);

% The columns represents the SNR (mu) value
Xrw = zeros(L+1, 10);
Xou = zeros(L+1, 10);
Xem = zeros(L+1, 10);

% BM w/ drift    
for mu = 1:1:10
    Xrw(1:L+1,mu) = mu * t + B.get_path;
end

% Initialize start value for ZRP particle in random medium
Xem(1, 1:10) = zeros(1, 10);

% Initialize start value for the OU process
Xou(1, 1:10) = ones(1, 10);

% Fill in
for j = 1:L
    dBt = sqrt(dt) * randn;
    
    for mu = 1:1:10
        spdchg = psi(j * dt, Xem(j, mu), 10 * mu);
    
        Xem(j+1, mu) = Xem(j,mu) + sqrt(spdchg) * dBt ... 
            - W.d_dx(Xem(j, mu), epsilon) * spdchg * dt;

        Xou(j+1,mu) = Xou(j,mu) + (1/(mu + 2)) * Xou(j,mu) * dt + dBt; 
    end
end

if (ISWRITE)
    % Write the medium into a CSV file
    str = sprintf('W_seed_%d.csv', seednum);
    writematrix(W.get_path(), str);
    
    % Write the particle dynamics into a CSV file
    str = sprintf('Xrw_L_%d_T_%d_seed_%d.csv', L, T, seednum);
    writematrix(Xrw, str);
    str = sprintf('Xou_L_%d_T_%d_seed_%d.csv', L, T, seednum);
    writematrix(Xou, str);
    str = sprintf('Xem_L_%d_T_%d_seed_%d.csv', L, T, seednum);
    writematrix(Xem, str);
end

% Plot the BM drift signal
if (ISPLOT)
    figure;
    xlabel('t', FontSize=12)
    ylabel('Xrw', 'FontWeight', 'normal')
    legend('Location','best')
    grid on
    hold on
    for mu = 1:1:10
        plot([0:dt:T], Xrw(1:L+1,mu), '-', 'LineWidth', 1, ...
            'DisplayName', num2str(mu))
    end
    hold off
    
    % Plot the LPSX signal
    figure;
    xlabel('t', FontSize=12)
    ylabel('Xzrp', 'FontWeight', 'normal')
    legend('Location','best')
    grid on
    hold on
    for mu = 1:1:10
        plot([0:dt:T], Xem(1:L+1,mu), '-', 'LineWidth', 1, ...
            'DisplayName', num2str(mu * 10))
    end
    hold off
    
    % Plot the OU signal
    figure;
    xlabel('t', FontSize=12)
    ylabel('Xou L', 'FontWeight', 'normal')
    legend('Location','best')
    grid on
    hold on
    for mu = 1:1:2
        plot([0:dt:T], Xou(1:L+1,mu), '-', 'LineWidth', 1, ...
            'DisplayName', num2str(1/ (mu + 2)))
    end
    hold off
    figure;
    xlabel('t', FontSize=12)
    ylabel('Xou M', 'FontWeight', 'normal')
    legend('Location','best')
    grid on
    hold on
    for mu = 3:1:5
        plot([0:dt:T], Xou(1:L+1,mu), '-', 'LineWidth', 1, ...
            'DisplayName', num2str(1/ (mu + 2)))
    end
    hold off
    figure;
    xlabel('t', FontSize=12)
    ylabel('Xou S', 'FontWeight', 'normal')
    legend('Location','best')
    grid on
    hold on
    for mu = 6:1:10
        plot([0:dt:T], Xou(1:L+1,mu), '-', 'LineWidth', 1, ...
            'DisplayName', num2str(1/ (mu + 2)))
    end
    hold off
end