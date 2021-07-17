function eps_r = eps_silicon(omega)
%EPS_SILICON Summary of this function goes here
%   Detailed explanation goes here
% https://arxiv.org/pdf/1705.05218.pdf
eps_inf = 0.81568;
hbar = 6.582119569e-16; % eV*s
sigma = [1.6934+2.084j 5.2573+8.0106j -1.7164+5.9939j -0.00528+0.32911j -3.8438+6.9298j];
Omega = [3.3736-0.11402j 3.6519-0.52378j 4.2877-0.21116j 5.3188-0.18434j 5.5064-1.7892j];
eps_r = eps_inf;
for i=1:5
    eps_r = eps_r + 1j.*sigma(i)./(omega.*hbar - Omega(i)) + 1j.*conj(sigma(i))./(omega.*hbar + conj(Omega(i)));
end
end

