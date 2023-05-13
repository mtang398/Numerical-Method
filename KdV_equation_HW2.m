L = 60;
c = 1;

N = 999;
h = 2.*L ./ N; % grid spacing
x = (0:N-1)*h - L; % include left endpoint

phi_sol = c/2 * sech((sqrt(c)/2)*(x)).^2; %actual Phi

phi_sol_grid = gridConstruct(81, L);
phi_sol_interp = MyInterpftPadded(phi_sol_grid, 81, N);

plot(x,phi_sol)
hold on
plot(x,phi_sol_interp)
legend('actual \phi_{sol}', 'interpreted \phi_{sol}')
title('Fourier Interpolation of \phi_{sol}')

%% For phi_sol function's error
modes = 27:20:487;
L_1_Error = [];
L_2_Error = [];
L_inf_Error = [];
for k = modes
    phi_sol_grid = gridConstruct(k, L);
    phi_sol_interp = MyInterpftPadded(phi_sol_grid, k, N);
    L_1_Error(end+1) = function_1_norm((phi_sol - phi_sol_interp.'), N, h) ./ function_1_norm(phi_sol, N, h);
    L_2_Error(end+1) = function_2_norm((phi_sol - phi_sol_interp.'), N, h) ./ function_2_norm(phi_sol, N, h);
    L_inf_Error(end+1) = max(abs(phi_sol - phi_sol_interp.')) ./ max(abs(phi_sol));
end

semilogy(modes, L_1_Error, '-o')
hold on
semilogy(modes, L_2_Error, '-o')
hold on
plot(modes, L_inf_Error, '-o')
legend('L_1 function norm', 'L_2 function norm', 'L_{\infty} function norm')
title('Error Function Norm Comparison')

%% For triangle wave
P = 10;
h = 2.*P ./ N; % grid spacing
x = (0:N-1)*h - P; % include left endpoint
tri = ((10+x).*(x<0 & x>=-10) + (10-x).*(x>=0 & x<=10)); % actual tri

tri_grid = triGridConstruct(81,P);
tri_interp = MyInterpftPadded(tri_grid,81,N);

plot(x, tri)
hold on
plot(x, tri_interp)
legend('actual triangle', 'interpreted triangle')
title('Fourier Interpolation of triangle wave')

%%
% for triangle wave's Error
modes = 27:20:487;
L_1_Error = [];
L_2_Error = [];
L_inf_Error = [];
for k = modes
    tri_grid = triGridConstruct(k, P);
    tri_interp = MyInterpftPadded(tri_grid, k, N);
    L_1_Error(end+1) = function_1_norm((tri - tri_interp.'), N, h) ./ function_1_norm(tri, N, h);
    L_2_Error(end+1) = function_2_norm((tri - tri_interp.'), N, h) ./ function_2_norm(tri, N, h);
    L_inf_Error(end+1) = max(abs(tri - tri_interp.')) ./ max(abs(tri));
end

semilogy(modes, L_1_Error, '-o')
hold on
semilogy(modes, L_2_Error, '-o')
hold on
plot(modes, L_inf_Error, '-o')
legend('L_1 function norm', 'L_2 function norm', 'L_{\infty} function norm')
title('Triangle Wave Error Function Norm')

%% For sawtooth wave
sawtooth = (10 + x) .* (x < 0 & x >= -10);  % actual sawtooth

sawtooth_grid = sawtoothGridConstruct(81,P);
sawtooth_interp = MyInterpftPadded(sawtooth_grid,81,N);

plot(x, sawtooth)
hold on
plot(x, sawtooth_interp)
legend('actual triangle', 'interpreted triangle')
title('Fourier Interpolation of sawtooth wave')

%% For sawtooth wave's error
modes = 27:20:487;
L_1_Error = [];
L_2_Error = [];
L_inf_Error = [];
for k = modes
    sawtooth_grid = sawtoothGridConstruct(k, P);
    sawtooth_interp = MyInterpftPadded(sawtooth_grid, k, N);
    L_1_Error(end+1) = function_1_norm((sawtooth - sawtooth_interp.'), N, h) ./ function_1_norm(sawtooth, N, h);
    L_2_Error(end+1) = function_2_norm((sawtooth - sawtooth_interp.'), N, h) ./ function_2_norm(sawtooth, N, h);
    L_inf_Error(end+1) = max(abs(sawtooth - sawtooth_interp.')) ./ max(abs(sawtooth));
end

semilogy(modes, L_1_Error, '-o')
hold on
semilogy(modes, L_2_Error, '-o')
hold on
plot(modes, L_inf_Error, '-o')
legend('L_1 function norm', 'L_2 function norm', 'L_{\infty} function norm')
title('Sawtooth Wave Error Function Norm')

%% Triangle Fourier Spectrum
f = triGridConstruct(163, P);
tri_hat = fftshift(fft(f)); % normal ordering
k = -81:81;
f_hat = abs(tri_hat);

semilogy(k, f_hat)
title('Triangle Wave Spectrum')

%% Sawtooth Fourier Spectrum
g = sawtoothGridConstruct(163, P);
sawtooth_hat = fftshift(fft(g)); % normal ordering
g_hat = abs(sawtooth_hat);

semilogy(k, g_hat);
title('Sawtooth Wave Spectrum')

%% Evaluate r.h.s of KdV
h = 2.*L ./ N; % grid spacing
x = (0:N-1)*h - L; % include left endpoint
rhs_KdV = 0.5 .* (c.^(5./2).*sinh(0.5.*sqrt(c)*x))./((cosh(0.5 .*sqrt(c).*x)).^3);

plot(x,rhs_KdV)
grid on
hold on

F = KdV_Spectral(N, 80,L);
plot(x,F)
legend('actual rhs', 'approximated rhs')

%% r.h.s Error

modes = 20:20:560;
L_1_Error = [];
L_2_Error = [];
L_inf_Error = [];
for k = modes
    kdv_interp = KdV_Spectral(N, k, L);
    L_1_Error(end+1) = function_1_norm((rhs_KdV - kdv_interp.'), N, h) ./ function_1_norm(rhs_KdV, N, h);
    L_2_Error(end+1) = function_2_norm((rhs_KdV - kdv_interp.'), N, h) ./ function_2_norm(rhs_KdV, N, h);
    L_inf_Error(end+1) = max(abs(rhs_KdV - kdv_interp.')) ./ max(abs(rhs_KdV));
end

semilogy(modes, L_1_Error, '-o')
hold on
semilogy(modes, L_2_Error, '-o')
hold on
plot(modes, L_inf_Error, '-o')
legend('L_1 function norm', 'L_2 function norm', 'L_{\infty} function norm')
title('r.h.s KdV Error Function Norm')

%% digits computation
M = 300:2:400;
new_L_2_Error = [];
for k = M
    kdv_interp = KdV_Spectral(N, k, L);
    new_L_2_Error(end+1) = function_2_norm((rhs_KdV - kdv_interp.'), N, h);
end

semilogy(M, new_L_2_Error, '-o')


function F = KdV_Spectral(N, n, L)
    phi = gridConstruct(n,L);
    h1 = 2.*L ./ n; % grid spacing
    x1 = (0:n-1)*h1 - L; % include left endpoint

    k = [0:n/2-1, 0, -n/2+1:-1] .* ( pi ./ L);
    phi_hat = (fft(phi));
    phi2_hat = (fft(phi.^2));
    rhs = ( 1i.*k.^3.*(phi_hat) - 3.*1i.*k.*fft(real(ifft(phi2_hat))));
    rhs_padd = [rhs(1:(n)/2).'; zeros(N-n,1); rhs((n)/2 + 1: n).'];
    F = real(ifft(rhs_padd)) .* (N/n);
end

function E1 = function_1_norm(f,N,h)
    E = 0;
    for c=1:N
        E = E + abs(f(c));
    end
    E1 = h .* E;
end

function E2 = function_2_norm(f,N,h)
    E = 0;
    for c=1:N
        E = E + (f(c)).^2;
    end
    E2 = sqrt(h .* E);
end

function psi_interp = MyInterpftPadded(psi, n, N)
    psi_hat = fft(psi); % get fourier coefficient for psi
    psi_padd = [psi_hat(1:(n+1)/2).'; zeros(N-n,1); psi_hat((n+1)/2 + 1: n).']; % padd zeros
    psi_interp = real(ifft((psi_padd)))* (N/n); % use IFFT to transform back to the real space
end

function psi = gridConstruct(n, L)
    h = 2.*L ./ n; % grid spacing
    x = (0:n-1)*h - L; % include left endpoint
    psi = 1/2 * sech((sqrt(1)/2)*(x)).^2;
end

function tri = triGridConstruct(n, P)
    h = 2.*P ./ n; % grid spacing
    x = (0:n-1)*h - P; % include left endpoint
    tri = ((10+x).*(x<0 & x>=-10) + (10-x).*(x>=0 & x<=10));
end

function sawtooth = sawtoothGridConstruct(n, P)
    h = 2.*P ./ n; % grid spacing
    x = (0:n-1)*h - P; % include left endpoint
    sawtooth = (10 + x) .* (x < 0 & x >= -10);
end