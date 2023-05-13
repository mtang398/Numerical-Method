%% 1.1(a)

% Graph fourier interpolation with n = 9
N = 999;
h = 2.*pi ./ N; % grid spacing
x= (0:N-1)*h - pi; % include left endpoint

a = 3;
f = exp(a .* cos(x)) ./ (2 .* pi .* besseli(0,a)); % actual f

psi = gridConstruct(9);
psi_interp = MyInterpftPadded(psi,9,N);

plot(x,f)
hold on
plot(x,psi_interp)
title('Figure 1: Fourier Interpolation of f with n = 9')
xlabel('x')
ylabel('y')
legend('Actual f', 'Fourier Interpolation')

%% Graph error function eps(x) = |f(x) - phi(x)| when n = 27

psi1 = gridConstruct(27);
psi1_interp = MyInterpftPadded(psi1,27,N);
eps = abs(f - psi1_interp.');

plot(x,eps)
title('Figure 2: Error Function eps(x) = |f(x) - phi(x)|')
xlabel('x')
ylabel('y')

%% 1.1(b)

f2 = gridConstruct(55);
psi_hat = fftshift(fft(f2)); % normal ordering
k = -27:27;
f_hat = abs(psi_hat);

semilogy(k, f_hat)
title('Figure 3: |Fourier Spectrum| vs Node Index (normal ordering)')
xlabel('Fourier Node Index')
ylabel('Magnitude of Fourier Spectrum')

%% 1.2(a)
E2 = [];
for k = 1:6
    E2(end+1) = EuclideanError(f,k,N);
end
node = 1:6;
pts = 3 .^ node;

semilogy(pts,E2,'-o')
title('Figure 4: Euclidean Error E2 when n = 3^k')
xlabel('Node #')
ylabel('E2[\phi(x)]')

%% 1.3
df = -a .* sin(x) .* exp(a .* cos(x)) ./ (2 .* pi .* besseli(0,a)); % actual f'

dphi = spectralDerivative(N,8);

plot(x,df)
hold on
plot(x,dphi)
title('Figure 5: Spectral Derivatives when n = 8')
xlabel('x')
ylabel('y')
legend('Actual df', 'dphi')

%%
E2_dphi = [];
for k = 2:9
    E2_dphi(end+1) = EuclideanError_dphi(df,k,N);
end
node1 = 2:9;
pts1 = 2 .^ node1;

semilogy(pts1,E2_dphi,'-o')
title('Figure 6: Euclidean Error E2_dphi when n = 2^k')
xlabel('Node #')
ylabel('E2[dphi(x)]')

%% 1.4
func = @(x) exp(a .* cos(x)) ./ (2 .* pi .* besseli(0,a));
F = [];
for i = x
    F(end+1) = integral(func,0,i);
end

Phi = IntegralEval(N,8);

plot(x,F)
hold on
plot(x,Phi)
title('Figure 7: Fourier Approximated Integral')
xlabel('x')
ylabel('y')
legend('built-in function', '\int_0^x \phi(t) dt')

%%

ActualPhi = IntegralEval(N,512);
Phi = IntegralEval(N,8);

plot(x,ActualPhi)
hold on
plot(x,Phi)
title('Figure 8: Fourier Approximated Integral')
xlabel('x')
ylabel('y')
legend('Approximated \int_0^x f(t) dt with large n', '\int_0^x \phi(t) dt')

%%
E2_F = []
for k = 2:8
    E2_F(end+1) = EuclideanError_dphi(df,k,N);
end
node2 = 2:8;
pts2 = 2 .^ node2;

semilogy(pts2,E2_F,'-o')
title('Figure 9: Euclidean Error E2[F] when n = 2^k')
xlabel('Node #')
ylabel('E2[F(x)]')

function E2_F = EuclideanError_F(ActualPhi, k, N)
    E = 0;
    m = 2.^k;
    h = 2.*pi ./ N;
    Phi = IntegralEval(N,m);
    
    for c=1:N
        E = E + (ActualPhi(c) - Phi(c)).^2;
    end
    E2_F = sqrt(h .* E);
end

function Phi = IntegralEval(N,n)
    h = 2.*pi ./ n; % grid spacing
    x = (0:n-1)*h - pi; % include left endpoint
    phi = exp(3 .* cos(x)) ./ (2 .* pi .* besseli(0,3));
    f0 = (fft(phi));
    ik = [0 -1i./(1:n/2-1) 0 -1i./(-n/2+1:-1)];
    phi_hat = f0 .* ik;

    H = 2.*pi ./ N; % grid spacing
    X = (0:N-1)*H - pi; % include left endpoint
    Phi_padd = [phi_hat(1:(n)/2).'; zeros(N-n,1); phi_hat((n)/2 + 1: n).'];
    Phi = (f0(1) .*X /n).'+ (ifft(Phi_padd)-sum(Phi_padd)) .* (N ./ n);
end

function E2_dphi = EuclideanError_dphi(df, k, N)
    E = 0;
    m = 2.^k;
    h = 2.*pi ./ N;
    dphi = spectralDerivative(N,m);
    
    for c=1:N
        E = E + (df(c) - dphi(c)).^2;
    end
    E2_dphi = sqrt(h .* E);
end

function dphi = spectralDerivative(N,n)
    phi = gridConstruct(n);
    %dphi = -3 .* sin(x) .* phi;
    phi_hat = fft(phi);
    k = [0:(n/2-1) 0 (-n/2 + 1):-1];
    ik = 1i .* k;
    dphi_hat = ik .* phi_hat;
    dphi_padd = [dphi_hat(1:(n)/2).'; zeros(N-n,1); dphi_hat((n)/2 + 1: n).'];
    dphi = real(ifft(dphi_padd))* (N/n);
end

function E2 = EuclideanError(f, k, N)
    E = 0;
    m = 3.^k;
    h = 2.*pi ./ N;
    psi = gridConstruct(m);
    psi_interp = MyInterpftPadded(psi,m,N);
    
    for c=1:N
        E = E + (f(c) - psi_interp(c)).^2;
    end
    E2 = sqrt(h .* E);
end

function psi_interp = MyInterpftPadded(psi, n, N)
    psi_hat = fft(psi); % get fourier coefficient for psi
    psi_padd = [psi_hat(1:(n+1)/2).'; zeros(N-n,1); psi_hat((n+1)/2 + 1: n).']; % padd zeros
    psi_interp = real(ifft((psi_padd)))* (N/n); % use IFFT to transform back to the real space
end

function psi = gridConstruct(n)
    h = 2.*pi ./ n; % grid spacing
    x = (0:n-1)*h - pi; % include left endpoint
    psi = exp(3 .* cos(x)) ./ (2 .* pi .* besseli(0,3));
end
