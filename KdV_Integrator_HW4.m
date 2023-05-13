% Parameters
L = 60;
c = 1;

% Actual phi
N = 1000;
h = L/N; % grid spacing
x = (0:N-1)*h - L/2; % include left endpoint

phi_sol = c/2 * sech((sqrt(c)/2)*(x)).^2; %actual Phi

dt = 0.01;
t0 = 0;
tf = L/c;

[t, phi] = SBDF2(N, 32, dt, t0, tf, tf/dt+1);
[t1, phi1] = SBDF2(N, 32, dt, t0, tf, 1);
[t2, phi2] = ETDRK2(N, 32, dt, t0, tf, tf/dt+1);
[t3, phi3] = ETDRK2(N, 32, dt, t0, tf, 1);

figure;
subplot(2,2,1);
plot(x, phi1, 'b-', 'LineWidth', 1.5)
hold on
plot(x, phi_sol, "--", 'MarkerSize', 10)
title('SBDF2 at t = 0')
xlabel('x')
ylabel('phi and phi_{sol}')
legend('phi', 'phi_{sol}', 'Location', 'best')
grid on

subplot(2,2,2);
plot(x, phi3, 'b-', 'LineWidth', 1.5)
hold on
plot(x, phi_sol, "--", 'MarkerSize', 10)
title('ETDRK2 at t = 0')
xlabel('x')
ylabel('phi and phi_{sol}')
legend('phi', 'phi_{sol}', 'Location', 'best')
grid on

subplot(2,2,3);
plot(x, phi, 'b-', 'LineWidth', 1.5)
hold on
plot(x, phi_sol, "--", 'MarkerSize', 10)
title(['SBDF2 at t = ' num2str(tf)])
xlabel('x')
ylabel('phi and phi_{sol}')
legend('phi', 'phi_{sol}', 'Location', 'best')
grid on

subplot(2,2,4);
plot(x, phi2, 'b-', 'LineWidth', 1.5)
hold on
plot(x, phi_sol, "--", 'MarkerSize', 10)
title(['ETDRK2 at t = ' num2str(tf)])
xlabel('x')
ylabel('phi and phi_{sol}')
legend('phi', 'phi_{sol}', 'Location', 'best')
grid on

sgtitle(['Comparison of SBDF2 and ETDRK2 methods at t=0 and t=' num2str(tf)])

%% 1.2.1
% Parameters
L = 60;
c = 1;
t0 = 0;
tf = 4*L/c;

% Grid sizes (number of Fourier modes)
grid_sizes = [32, 64, 128, 256];

% Time step sizes
dts = [0.01, 0.02, 0.04, 0.08];
% Create separate figures for SBDF2 and ETDRK2
figure;
sgtitle('SBDF2 method for different grid sizes and time step sizes');
figure;
sgtitle('ETDRK2 method for different grid sizes and time step sizes');

% Loop over time step sizes
for dt_idx = 1:length(dts)
    dt = dts(dt_idx);
    T = L/dt+1;
    % Plot phi_sol
    x = (0:N-1)*h - L/2;
    phi_sol = c/2 * sech((sqrt(c)/2)*(x)).^2;

    % Create subplots for each method and time step size
    figure(1);
    subplot(3, 2, dt_idx);
    hold on;
    plot(x, phi_sol, 'k--', 'DisplayName', 'phi_{sol}');
    
    figure(2);
    subplot(3, 2, dt_idx);
    hold on;
    plot(x, phi_sol, 'k--', 'DisplayName', 'phi_{sol}');
    
    % Loop over grid sizes
    for grid_idx = 1:length(grid_sizes)
        n = grid_sizes(grid_idx);
        
        % Simulate the KdV equation
        [~, phi_SBDF2] = SBDF2(N, n, dt, t0, tf, T);
        [~, phi_ETDRK2] = ETDRK2(N, n, dt, t0, tf, T);
        
        % Plot the SBDF2 results
        figure(1);
        plot(x, phi_SBDF2, 'DisplayName', sprintf('SBDF2 (n = %d)', n));
        xlabel('x');
        ylabel('phi');
        title(sprintf('dt = %.2f', dt));
        if dt_idx == 1
            legend('Location','best')
        end
        grid on;
        
        % Plot the ETDRK2 results
        figure(2);
        plot(x, phi_ETDRK2, 'DisplayName', sprintf('ETDRK2 (n = %d)', n));
        xlabel('x');
        ylabel('phi');
        title(sprintf('dt = %.2f', dt));
        if dt_idx == 1
            legend('Location','best')
        end
        grid on;
    end
    
    % Release hold on subplots
    figure(1);
    hold off;
    
    figure(2);
    hold off;
end


%% 1.2.2
% Parameters
L = 60;
c = 1;

% Actual phi
N = 1000;
n = 128;
h = L/N; % grid spacing
x = (0:N-1)*h - L/2; % include left endpoint

% Time step sizes
k = -2:5;
dt_values = 0.01 ./ 2.^k;

% Initialize error storage
error_SBDF2 = zeros(length(dt_values), 1);
error_ETDRK2 = zeros(length(dt_values), 1);
slope_SBDF2 = zeros(length(dt_values)-1, 1);
slope_ETDRK2 = zeros(length(dt_values)-1, 1);

% Loop over time step sizes
for i = 1:8
    dt = dt_values(i);
    dt_half = dt/2;
    
    % SBDF2 method
    [~, phi_dt] = SBDF2(N, n, dt, t0, 15, 15/dt+1);
    [~, phi_dt_half] = SBDF2(N, n, dt_half, t0, 15, 15/dt_half+1);
    
    % Error calculation for SBDF2
    error_SBDF2(i) = norm(phi_dt - phi_dt_half, "inf");
    
    % ETDRK2 method
    [~, phi_dt_ETDRK2] = ETDRK2(N, n, dt, t0, 15, 15/dt+1);
    [~, phi_dt_half_ETDRK2] = ETDRK2(N, n, dt_half, t0, 15, 15/dt_half+1);
    
    % Error calculation for ETDRK2
    error_ETDRK2(i) = norm(phi_dt_ETDRK2 - phi_dt_half_ETDRK2, "inf");
end

% Plot errors and second-order convergence
figure;
loglog(dt_values, error_SBDF2, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
loglog(dt_values, error_ETDRK2, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10);

% Second-order convergence line
loglog(dt_values, dt_values.^2, 'k--', 'LineWidth', 1.5);

% Plot settings
xlabel('Time step size (\Delta t)');
ylabel('Error');
title('Error vs. Time step size for SBDF2 and ETDRK2');
legend('SBDF2', 'ETDRK2', 'Second-order convergence', 'Location', 'best');
grid on;

% Compute slopes
for i = 1:length(dt_values)-1
    slope_SBDF2(i) = diff(log(error_SBDF2(i:i+1))) / diff(log(dt_values(i:i+1)));
    slope_ETDRK2(i) = diff(log(error_ETDRK2(i:i+1))) / diff(log(dt_values(i:i+1)));
end

% Output slopes
disp('Slopes:');
disp(['SBDF2: ' num2str(slope_SBDF2')]);
disp(['ETDRK2: ' num2str(slope_ETDRK2')]);

%% 1.2.3
% Parameters
L = 60;
c = 1;

% Actual phi
N = 1000;
n = 128;
h = L/N; % grid spacing
x = (0:N-1)*h - L/2; % include left endpoint

% Time step sizes
k = -2:5;
dt_values = 0.01 ./ 2.^k;

% Initialize error storage
error_SBDF2 = zeros(length(dt_values), 1);
error_ETDRK2 = zeros(length(dt_values), 1);
slope_SBDF2 = zeros(length(dt_values)-1, 1);
slope_ETDRK2 = zeros(length(dt_values)-1, 1);

% Loop over time step sizes
for i = 1:8
    dt = dt_values(i);
    
    % SBDF2 method
    [~, phi_dt] = SBDF2(N, n, dt, t0, L, L/dt+1);
    
    % Error calculation for SBDF2
    error_SBDF2(i) = norm(phi_dt - phi_sol, 1);
    
    % ETDRK2 method
    [~, phi_dt_ETDRK2] = ETDRK2(N, n, dt, t0, L, L/dt+1);
    
    % Error calculation for ETDRK2
    error_ETDRK2(i) = norm(phi_dt_ETDRK2 - phi_sol, 1);
end

% Plot errors and second-order convergence
figure;
loglog(dt_values, error_SBDF2, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
loglog(dt_values, error_ETDRK2, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 10);

% Second-order convergence line
loglog(dt_values, dt_values.^2, 'k--', 'LineWidth', 1.5);

% Plot settings
xlabel('Time step size (\Delta t)');
ylabel('Error');
title('Error vs. Time step size for SBDF2 and ETDRK2');
legend('SBDF2', 'ETDRK2', 'Second-order convergence', 'Location', 'best');
grid on;

% Compute slopes
for i = 1:length(dt_values)-1
    slope_SBDF2(i) = diff(log(error_SBDF2(i:i+1))) / diff(log(dt_values(i:i+1)));
    slope_ETDRK2(i) = diff(log(error_ETDRK2(i:i+1))) / diff(log(dt_values(i:i+1)));
end

% Output slopes
disp('Slopes:');
disp(['SBDF2: ' num2str(slope_SBDF2')]);
disp(['ETDRK2: ' num2str(slope_ETDRK2')]);

%% 1.2.4 
% Parameters
L = 60;
c = 1;
t0 = 0;
tf = 4*L/c;

% Grid sizes (number of Fourier modes)
grid_sizes = [32, 64, 128];

% Time step sizes
dts = 0.04;
% Create separate figures for SBDF2 and ETDRK2
figure;
sgtitle('SBDF2, dt = 0.04');
figure;
sgtitle('ETDRK2, dt = 0.04');

% Loop over time step sizes
for dt_idx = 1:length(dts)
    dt = dts(dt_idx);
    T = L/dt+1;
    % Plot phi_sol
    x = (0:N-1)*h - L/2;
    phi_sol = c/2 * sech((sqrt(c)/2)*(x)).^2;

    % Create subplots for each method and time step size
    figure(1);
    hold on;
    plot(x, phi_sol, 'k--', 'DisplayName', 'phi_{sol}');
    
    figure(2);
    hold on;
    plot(x, phi_sol, 'k--', 'DisplayName', 'phi_{sol}');
    
    % Loop over grid sizes
    for grid_idx = 1:length(grid_sizes)
        n = grid_sizes(grid_idx);
        
        % Simulate the KdV equation
        [~, phi_SBDF2] = SBDF2(N, n, dt, t0, tf, T);
        [~, phi_ETDRK2] = ETDRK2(N, n, dt, t0, tf, T);
        
        % Plot the SBDF2 results
        figure(1);
        plot(x, phi_SBDF2, 'DisplayName', sprintf('SBDF2 (n = %d)', n));
        xlabel('x');
        ylabel('phi');
        if dt_idx == 1
            legend('Location','best')
        end
        grid on;
        
        % Plot the ETDRK2 results
        figure(2);
        plot(x, phi_ETDRK2, 'DisplayName', sprintf('ETDRK2 (n = %d)', n));
        xlabel('x');
        ylabel('phi');
        if dt_idx == 1
            legend('Location','best')
        end
        grid on;
    end
    
    % Release hold on subplots
    figure(1);
    hold off;
    
    figure(2);
    hold off;
end



function c_n = coeff(n)
    L = 60; 
    h = L/n; % grid spacing
    nodes = (0:n-1)*h - L/2; % include left endpoint
    c_n = fft(phi_solKDV(nodes, 0));
end

function phi = phi_solKDV(x, t)
    c = 1; 
    phi = c/2 * sech((sqrt(c)/2)*(x - c*t)).^2;
end

function c_n_padded = coeff_padd(N, n)
    padd_n = N - n;
    c_n = coeff(n);
    half_length = floor(length(c_n) / 2);
    c_n_padded = [c_n(1:half_length), zeros(1, padd_n), c_n(half_length+1:end)];
end

function A = Construct_A(N)
    L = 60;
    k_big = (2 * pi / L) * fftshift((-floor(N/2):floor(N/2)-1));
    A = 1i * (k_big.^3);
end

function B_result = B(N, phihat, n)
    L = 60;
    k_big = (2 * pi / L) * fftshift((-floor(N/2):floor(N/2)-1));
    phisquarehat = (N / n) * fft(real(ifft(phihat).^2));
    B_result = -3 * 1i * k_big .* phisquarehat;
end

function [t, phi] = SBDF2(N, n, dt, t0, tf, T)
    % Initialize the time vector and phi matrix
    t = t0:dt:tf;
    % k = ((N/n) * pi / L) * fftshift((-floor(n/2):floor(n/2)-1));
    A = Construct_A(N);
    phi_hat_0 = coeff_padd(N,n);
    phihat = zeros(length(t), length(phi_hat_0), 'like', complex(0));
    phihat(1, :) = phi_hat_0;
    phihat(2, :) = phihat(1, :) + dt * (A .* phihat(1, :) + B(N, phihat(1, :), n));

    A_const = 1 ./ (ones(size(A)) - (2*dt/3) .* A);

    for i = 2:length(t) - 1
        phihat(i+1, :) = A_const .* ((4/3)*phihat(i, :) - (1/3) * phihat(i-1, :) ...
            + (2*dt/3)*(2*B(N, phihat(i,:), n) - B(N, phihat(i-1,:), n)));
    end

    phi = (N/n) * real(ifft(phihat(T,:)));
end

function [t, phi] = ETDRK2(N, n, dt, t0, tf, T)
    % Initialize the time vector and phi matrix
    t = t0:dt:tf;
    A = Construct_A(N);
    A1 = A(A ~= 0);
    A_inv = 1 ./ A1;
    A_inv = [0; A_inv(:)].';
    A_inv2 = 1 ./ (A1 .^ 2);
    A_inv2 = [0; A_inv2(:)].';
    A_exp = exp(A .* dt);

    phi_hat_0 = coeff_padd(N, n);
    phihat = zeros(length(t), length(phi_hat_0), 'like', complex(0));
    phihat(1, :) = phi_hat_0;

    predconst = A_inv .* (A_exp - ones(size(A_exp)));
    corrconst = A_inv2 .* (A_exp - ones(size(A_exp)) - A .* dt) / dt;

    % Time-stepping loop
    for i = 1:length(t) - 1
        phihatn = phihat(i, :);
        Bphihatn = B(N, phihatn, n);
        predictor = A_exp .* phihatn + predconst .* Bphihatn;
        phihat(i+1, :) = predictor + corrconst.* (B(N, predictor, n) - Bphihatn);
    end

    phi = (N/n) * real(ifft(phihat(T,:)));
end