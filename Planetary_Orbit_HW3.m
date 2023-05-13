% Constant
G = 1;

% Parameters
r = 1;  % radius of orbit
w = 1;  % angular frequency of orbit
r0 = [r.*cos(0), r.*sin(0)];
v0 = [0, w.*r0(1)];
T = 2.*pi./w;
m0 = r.^3 .* w.^2 ./ G;

% Number of time steps
M = 1000;
dt = T/M; % time step size

[p, v, t] = integrate_circular_orbit(dt, M, r0, v0, m0, G);

% Plot the trajectory
figure;
plot(p(1,:), p(2,:), 'LineWidth', 1.5);
hold on;
plot(0, 0, 'ro', 'MarkerSize', 8); % Position of the star
title('Circular Orbit');
xlabel('x');
ylabel('y');
axis equal;
grid on;
hold off;

%%
% Error analysis
k = [5,6,7,8,9,10,11,12,13,14,15];
M_values = 2.^k;
errors = zeros(1, length(M_values));

for j = 1:length(M_values)
    M = M_values(j);
    dt = T / M;
    
    [p, v, t] = integrate_circular_orbit(dt, M, r0, v0, m0, G);
    
    % Compute the error
    error_max = Error_BS_RK23(r,w,p,v,t);
    errors(j) = error_max;
end

% Plot the error
figure;
loglog(T ./ M_values, errors, '-o', 'LineWidth', 1.5);
title('Error vs. Time Step Size');
xlabel('Time Step Size');
ylabel('Error');
grid on;

% Calculate the order of accuracy
order = diff(log(errors)) ./ diff(log(T ./ M_values));
fprintf('Order of accuracy: %f\n', mean(order));

%%

M_val_test = 10:32;
errors = [];
for k = M_val_test
    dt = T/k;

    [p, v, t] = integrate_circular_orbit(dt, k, r0, v0, m0, G);
    % Exact solution
    x_exact = r * cos(w * t);
    y_exact = r * sin(w * t);

    % Compute the error
    error_x = abs(x_exact - p(1, :));
    error_y = abs(y_exact - p(2, :));
    error_max = max(max(error_x), max(error_y));
    errors(end+1) = error_max;
    if error_max < 0.01
        fprintf('Smallest Possible Step Size Number: %f\n', k);
        fprintf('Largest Possible Time Step Size: %f\n', dt);
        break
    end
end

Z = 10:k;

% Plot the error
figure;
plot(T ./ Z, errors.', '-o', 'LineWidth', 1.5);
title('Error vs. Time Step Size 2');
xlabel('Time Step Size');
ylabel('Error');
grid on;

%%
% Solution 1 when M = 22
[p1, v1, t1] = integrate_circular_orbit(T/22, 22, r0, v0, m0, G);
% Solution 2 when M = 44
[p2, v2, t2] = integrate_circular_orbit(T/44, 44, r0, v0, m0, G);
% Richardson Extrapolation
% Parameters
M_1 = 22;
M_2 = 2 * M_1;
dt_1 = T / M_1;
dt_2 = T / M_2;

% Integrate the orbit with h and h/2 time step sizes
[p_1, ~, t_1] = integrate_circular_orbit(dt_1, M_1, r0, v0, m0, G);
[p_2, ~, t_2] = integrate_circular_orbit(dt_2, M_2, r0, v0, m0, G);

% Richardson extrapolation
p_extrapolated = zeros(2, M_1+1);
for i = 1:(M_1+1)
    p_extrapolated(:, i) = (8 * p_2(:, 2*i-1) - p_1(:, i)) / 7;
end

% Plot the three solutions
figure;
hold on;
plot(p1(1,:), p1(2,:), 'LineWidth', 1.5, 'DisplayName', 'M = 22');
plot(p2(1,:), p2(2,:), 'LineWidth', 1.5, 'DisplayName', 'M = 44');
plot(p_extrapolated(1,:), p_extrapolated(2,:), 'LineWidth', 1.5, 'DisplayName', 'Richardson Extrapolated');
plot(0, 0, 'ro', 'MarkerSize', 8); % Position of the star
title('Circular Orbit - Comparison of Solutions');
xlabel('x');
ylabel('y');
axis equal;
grid on;
legend('show');
hold off;

%%
% Calculate the true solution
r_true1 = r * [cos(w * t1); sin(w * t1)];
r_true2 = r * [cos(w * t2); sin(w * t2)];
r_true_extrapolated = r * [cos(w * t_1); sin(w * t_1)];

% Calculate the errors
error1 = vecnorm(r_true1 - p1);
error2 = vecnorm(r_true2 - p2);
error_extrapolated = vecnorm(r_true_extrapolated - p_extrapolated);

% Plot the errors as a function of time
figure;
plot(t1, error1, '-o', 'LineWidth', 1.5);
hold on;
plot(t2, error2, '-x', 'LineWidth', 1.5);
plot(t_1, error_extrapolated, '-s', 'LineWidth', 1.5);
hold off;

title('Error vs. Time');
xlabel('Time (t)');
ylabel('Error (e^k)');
legend('Solution 1 (M = 22)', 'Solution 2 (M = 44)', 'Richardson Extrapolated Solution');
grid on;

% Plot the errors as a function of time
figure;
semilogy(t1, error1, '-o', 'LineWidth', 1.5);
hold on;
semilogy(t2, error2, '-x', 'LineWidth', 1.5);
semilogy(t_1, error_extrapolated, '-s', 'LineWidth', 1.5);
hold off;

title('Error vs. Time');
xlabel('Time (t)');
ylabel('Error (e^k)');
legend('Solution 1 (M = 22)', 'Solution 2 (M = 44)', 'Richardson Extrapolated Solution');
grid on;

%%
fprintf('Maximum Error of Richardson Extrapolation: %f\n', max(error_extrapolated));

%%
% Restating the parameters to make sure the code works
% Constant
G = 1;

% Parameters
r = 1;  % radius of orbit
w = 1;  % angular frequency of orbit
r0 = [r.*cos(0), r.*sin(0)];
v0 = [0, w.*r0(1)];
m0 = r.^3 .* w.^2 ./ G;

% Velocity scaling factors
factors = [1, 1/2, 1/4, 1/8];

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-2 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    %m0_scaled = r.^3 .* (w*factor).^2 ./ G;
    
    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);
    %[pos_exact, vel_exact] = exact_two_body_solution(r0, v0_scaled, m0, T_scaled, G);

    % Compute the error in the position after one orbit
    pos_err = norm(pos(:,1) - pos(:,end));
    fprintf('Error in position for velocity factor %.2f: %.2e\n', factor, pos_err);

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(pos(1,:), pos(2,:), '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('x position');
    ylabel('y position');
    title(sprintf('Orbits for Factor: %.2f', factor));
    counter = counter + 1;
end

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-2 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the error in the position after one orbit
    pos_err = norm(pos(:,1) - pos(:,end));
    fprintf('Error in position for velocity factor %.2f: %.2e\n', factor, pos_err);

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(time, pos(1,:), '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('time');
    ylabel('x position');
    title(sprintf('x(t) for Factor: %.2f', factor));
    counter = counter + 1;
end

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-2 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the number of time step
    fprintf('The number of time steps it take for factor %.2f: %.2f\n', factor, length(time));

    % Compute the error in the position after one orbit
    pos_err = norm(pos(:,1) - pos(:,end));
    fprintf('Error in position for velocity factor %.2f: %.2e\n', factor, pos_err);

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(time, dt_list, '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('time');
    ylabel('step size');
    title(sprintf('Step Size for Factor: %.2f', factor));
    counter = counter + 1;
end

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-2 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the number of time step
    fprintf('The number of time steps it take for factor %.2f: %.2f\n', factor, length(time));

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(time, error, '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('time');
    ylabel('step size');
    title(sprintf('Step Size for Factor: %.2f', factor));
    counter = counter + 1;
end
%%
% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-3 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    %m0_scaled = r.^3 .* (w*factor).^2 ./ G;
    
    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the error in the position after one orbit
    pos_err = norm(pos(:,1) - pos(:,end));
    fprintf('Error in position for velocity factor %.2f: %.2e\n', factor, pos_err);

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(pos(1,:), pos(2,:), '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('x position');
    ylabel('y position');
    title(sprintf('Orbits for Factor: %.2f', factor));
    counter = counter + 1;
end

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-3 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the error in the position after one orbit
    pos_err = norm(pos(:,1) - pos(:,end));
    fprintf('Error in position for velocity factor %.2f: %.2e\n', factor, pos_err);

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(time, pos(1,:), '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('time');
    ylabel('x position');
    title(sprintf('x(t) for Factor: %.2f', factor));
    counter = counter + 1;
end

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-3 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the error in the position after one orbit
    pos_err = norm(pos(:,1) - pos(:,end));
    fprintf('Error in position for velocity factor %.2f: %.2e\n', factor, pos_err);

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(time, dt_list, '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('time');
    ylabel('step size');
    title(sprintf('Step Size for Factor: %.2f', factor));
    counter = counter + 1;
end

% Main script to run the simulation and plot results
figure
hold on

% Loop over velocity scaling factors
counter = 1;
for factor = factors
    v0_scaled = v0 * factor;
    err_tol = 1e-3 * norm(r0);
    E = - (G*m0)/(norm(r0)) + norm(v0_scaled)^2 / 2;
    a = - (G*m0)/(2*E);
    T_scaled = 2 * pi * sqrt(a^3 / (G*m0));

    % Run the adaptive BS-RK3 integrator
    [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0_scaled, T_scaled, G, m0, err_tol);

    % Compute the number of time step
    fprintf('The number of time steps it take for factor %.2f: %.2f\n', factor, length(time));

    % Plot x(t) for all points along the integration
    subplot(2, 2, counter);
    plot(time, error, '-o','DisplayName', sprintf('Factor: %.2f', factor));
    xlabel('time');
    ylabel('step size');
    title(sprintf('Step Size for Factor: %.2f', factor));
    counter = counter + 1;
end
%% Animation

% Masses (in arbitrary units):
m1 = 1;
m2 = 1;
m3 = 1;
m4 = 3;
m5 = 2;

% Positions (in arbitrary units):
x1 = [0, 0];
x2 = [2, 0];
x3 = [0, 2];
x4 = [1, 1];
x5 = [-1, -1];

% Velocities (in arbitrary units):
v1 = [0, 0];
v2 = [0, 1];
v3 = [1, 0];
v4 = [0, -0.5];
v5 = [0.5, 0];

% Combine initial conditions into a single state vector
initial_conditions = [x1, x2, x3, x4, x5, v1, v2, v3, v4, v5]';

% Set up time span for integration
tspan = [0, 6]; % Time span for the simulation (in arbitrary time units)

h = 10 / 100;
tol = 1e-2*h / 10;
[t, y, h] = bogacki_shampine(@n_body, 0, 10, initial_conditions, h, tol);

% Set up the figure for movie creation
figure;
axis([-4, 4, -4, 4]);
axis equal;
hold on;

% Prepare the objects for the movie
body_colors = {'r', 'g', 'b', 'c', 'm'};
body_markers = {'o', 's', 'd', 'v', '^'};
body_handles = [];

for i = 1:5
    body_handles(i) = plot(y(1, 2*i-1), y(1, 2*i), body_markers{i}, 'Color', body_colors{i}, 'MarkerFaceColor', body_colors{i});
end

% Set up the movie structure
movie_frames = struct('cdata', cell(1, length(t)), 'colormap', cell(1, length(t)));
disp(length(t))
% Loop through the time steps and update the positions of the bodies
for k = 1:length(t)
    for i = 1:5
        set(body_handles(i), 'XData', y(k, 2*i-1), 'YData', y(k, 2*i));
    end
    drawnow;
    movie_frames(k) = getframe(gcf);
    %disp(k)
end

% Create the movie file (in AVI format)
movie_name = 'n_body_problem_simulation.avi';
movie_writer = VideoWriter(movie_name);
open(movie_writer);
writeVideo(movie_writer, movie_frames);
close(movie_writer);

% Define the function for the acceleration due to gravity
function a = gravity(p, m0, G)
    r = norm(p);
    a = -G * m0 / r^3 * p;
end

function [pos, vel, time, dt_list, error] = adaptive_bs_rk3(r0, v0, T, G, m0, tol)
    % Initialize arrays to store position, velocity, and time
    pos = [];
    vel = [];
    time = [];
    dt_list = [];
    error = [];
    
    % Initial values
    pos(:,1) = r0;
    vel(:,1) = v0;
    time(1) = 0;
    error(1) = 0;
    
    % Initialize step size
    dt = T / 22;
    dt_list(1) = dt;

    err_tol = tol*dt/T;
    
    while time(end) < T
        
        % Compute k1, k2, k3
        k1_v = dt * gravity(pos(:,end), m0, G);
        k1_r = dt * vel(:,end);
        
        k2_v = dt * gravity(pos(:,end) + 0.5*k1_r, m0, G);
        k2_r = dt * (vel(:,end) + 0.5*k1_v);
        
        k3_v = dt * gravity(pos(:,end) + 0.75*k2_r, m0, G);
        k3_r = dt * (vel(:,end) + 0.75*k2_v);
        
        % Calculate new position and velocity using BS-RK3 method
        r_new = pos(:,end) + (2*k1_r + 3*k2_r + 4*k3_r) / 9;
        v_new = vel(:,end) + (2*k1_v + 3*k2_v + 4*k3_v) / 9;
        
        % Compute k4
        k4_v = dt * gravity(r_new, m0, G);
        k4_r = dt * v_new;
        
        % Calculate the refined position and velocity
        r_refined = pos(:,end) + (7*k1_r + 6*k2_r + 8*k3_r + 3*k4_r) / 24;
        v_refined = vel(:,end) + (7*k1_v + 6*k2_v + 8*k3_v + 3*k4_v) / 24;
        
        % Estimate error
        pos_err = norm(r_refined - r_new, 'inf');
        vel_err = norm(v_refined - v_new, 'inf');
        err = max(pos_err, vel_err);
        
        % Check if error is within tolerance
        if err <= err_tol
            % Update position, velocity, and time
            dt_list(end+1) = dt;
            time(end+1) = time(end) + dt;
            pos(:,end+1) = r_new;
            vel(:,end+1) = v_new;
            error(end+1) = err;

            % Update step size for the next iteration
            dt = min(0.9 * dt * (err_tol / err)^(1/4), T - time(end));
        else
            % Update step size without updating time
            dt = min(0.9 * dt * (err_tol / err)^(1/3), T - time(end));
        end
    end
end


% Function to compute the maximum error of our approximation
function error_max = Error_BS_RK23(r ,w, p, v, t)
    % Exact solution
    x_exact = r * cos(w * t);
    y_exact = r * sin(w * t);

    % Compute the error
    error_x = abs(x_exact - p(1, :));
    error_y = abs(y_exact - p(2, :));
    error_max = max(sqrt(error_x.^2 + error_y.^2));
end

% Function to integrate the system using the BS-RK23 method
function [p, v, t] = integrate_circular_orbit(dt, M, r0, v0, m0, G)
    p = zeros(2, M+1);
    v = zeros(2, M+1);
    t = zeros(1, M+1);
    p(:,1) = r0;
    v(:,1) = v0;
    t(1) = 0;

    for i = 1:M
        t(i+1) = t(i) + dt;
        K1_v = gravity(p(:,i), m0, G);
        K1_p = v(:,i);
        K2_v = gravity(p(:,i) + 0.5*dt*K1_p, m0, G);
        K2_p = v(:,i) + 0.5*dt*K1_v;
        K3_v = gravity(p(:,i) + 0.75*dt*K2_p, m0, G);
        K3_p = v(:,i) + 0.75*dt*K2_v;
        p(:,i+1) = p(:,i) + dt*(2/9*K1_p + 1/3*K2_p + 4/9*K3_p);
        v(:,i+1) = v(:,i) + dt*(2/9*K1_v + 1/3*K2_v + 4/9*K3_v);
    end
end

function dydt = n_body(t, y)
    G = 1;
    m = [1, 1, 1, 3, 2];
    N = length(m);
    
    dydt = zeros(20, 1);
    
    for i = 1:N
        pos_i = y(2*i-1:2*i);
        vel_i = y(10+2*i-1:10+2*i);
        dydt(2*i-1:2*i) = vel_i;
        
        acceleration = [0; 0];
        for j = 1:N
            if i ~= j
                pos_j = y(2*j-1:2*j);
                r = pos_j - pos_i;
                acceleration = acceleration + G * m(j) * r / norm(r)^3;
            end
        end
        dydt(10+2*i-1:10+2*i) = acceleration;
    end
end

function [time, solution, stepSizes] = bogacki_shampine(func, tStart, tEnd, yStart, stepSize, tolerance)
    % Initialize arrays to store solution and time points
    time = tStart;
    solution = yStart;
    stepSizes = [];

    % Initialize current time and solution
    tCurr = tStart;
    yCurr = yStart;
    while tCurr < tEnd
        stepSizes = [stepSizes, stepSize];
        [est2, est3] = bs_step(func, tCurr, yCurr, stepSize);

        % Calculate the error estimate e
        errorEstimate = norm(est3 - est2, Inf);

        % Check if error is within tolerance
        if errorEstimate < tolerance
            stepSize = min(0.9 * stepSize * (tolerance / errorEstimate)^(1/4), (tEnd - tCurr));
            tCurr = tCurr + stepSize;
            yCurr = est2;
            time = [time, tCurr];
            solution = [solution, yCurr]; 
        end
        % If error is not within tolerance, reduce the step size
        stepSize = min(0.9 * stepSize * (tolerance / errorEstimate)^(1/3), (tEnd - tCurr));
    end

    % Convert arrays to column vectors
    time = time';
    solution = solution';
    stepSizes = stepSizes';

    % Bogacki-Shampine step function
    function [est2, est3] = bs_step(func, tCurr, yCurr, step)
        k1 = step * func(tCurr, yCurr);
        k2 = step * func(tCurr + step/2, yCurr + k1/2);
        k3 = step * func(tCurr + 3*step/4, yCurr + 3*k2/4);
        est2 = yCurr + (2*k1 + 3*k2 + 4*k3) / 9;
        k4 = step * func(tCurr + step, est2);
        est3 = yCurr + (7*k1/24 + k2/4 + k3/3 + k4/8);
    end
end

