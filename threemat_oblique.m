mu_0 = 1.256637061*10^-6; % H/m
eps_0 = 8.854187817*10^-12; % F/m
c = 299792458; % m/s
hbar = 6.582119569e-16; % eV*s
Z0 = sqrt(mu_0/eps_0);

N = 1000;
%lambda = linspace(400e-9, 1e-6, N); % m
%omega = (2*pi*c)./lambda; % 1/s
omega = 2*pi*linspace(5e12,30e12,N); % 1/s
lambda = (2*pi*c)./omega; % m
k0 = 2*pi./lambda; % 1/m

% delta_omega = omega(2) - omega(1);
% eps_2_re = 4.2 + zeros(size(omega));
% eps_2_im = zeros(size(omega));
% a=zeros(size(omega));
% b=zeros(size(omega));
% 
% alpha = 1;
% beta1=0;
% for k=2:N
% b(1)=beta1+eps_2_re(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(1)^2);
% beta1=b(1);
% end
% eps_2_im(1)=-2/pi*delta_omega*b(1)*omega(1)^(1-2*alpha);
% %First element of the output: the principal part integration
% %is computed by excluding the first element of the input
% 
% alpha1=0;
% for k=1:N-1
%     a(N)=alpha1+eps_2_re(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(N)^2);
%     alpha1=a(N);
% end
% eps_2_im(N)=-2/pi*delta_omega*a(N)*omega(N)^(1-2*alpha);
% %Last element of the output: the principal part integration
% %is computed by excluding the last element of the input.
% 
% for j=2:N-1
%     %Loop on the inner components of the output vector.
%     alpha1=0;
%     beta1=0;
%     for k=1:j-1
%         a(j)=alpha1+eps_2_re(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(j)^2);
%         alpha1=a(j);
%     end
%     for k=j+1:N
%         b(j)=beta1+eps_2_re(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(j)^2);
%         beta1=b(j);
%     end
%     eps_2_im(j)=-2/pi*delta_omega*(a(j)+b(j))*omega(j)^(1-2*alpha);
%     %Last element of the output: the principal part integration
%     %is computed by excluding the last element of the input
% end
% 
% figure;
% plot(omega/(2*pi*1e12),eps_2_re,'r'); hold on;
% plot(omega/(2*pi*1e12),eps_2_im,'b');

beta_1 = k0;
Z_1 = Z0;

N_theta = 1;
theta_incident = linspace(0, 0);%60/180*pi, N_theta);
D = 10e-6; % m
legend_title = cell(1,N_theta);
transmitted_TE = figure;
transmitted_TM = figure;
reflected_TE = figure;
reflected_TM = figure;

for i=1:N_theta
    theta_1 = theta_incident(i);
    Gamma_TE = zeros(size(omega));
    Gamma_TM = zeros(size(omega));
    T_TE = zeros(size(omega));
    T_TM = zeros(size(omega));
    legend_title{1,i} = "\theta_i = " + num2str(theta_1*180/pi) + "^o";
    for j=1:N
        o = omega(j);
        % get refractive index of silicon
        % using Drude-Lorentz model for Silicon from https://arxiv.org/pdf/1705.05218.pdf
        n_Si = sqrt(eps_silicon(o));
        beta_3 = k0(j).*3.5;%real(n_Si);
        Z_3 = Z0./3.5;

        % dielectric propagation constants
        eps_2 = 4.2*exp(1j*0.005);
        eps_2_re = real(eps_2);
        eps_2_im = imag(eps_2);
        alpha_2 = k0(j).*imag(sqrt(eps_2));
        beta_2 = k0(j).*real(sqrt(eps_2));
        Z_2 = Z0./sqrt(eps_2_re+1j*eps_2_im);

        cos_theta_2 = sqrt(beta_2.^2-beta_1(j).^2.*(sin(theta_1)).^2)./beta_2;
        cos_theta_3 = sqrt(beta_3.^2-beta_1(j).^2.*(sin(theta_1)).^2)./beta_3;
        TE_A = [1, -1, -1, 0;
             0, exp(-D/2.*alpha_2-1j.*D.*beta_2.*cos_theta_2), exp(D/2.*alpha_2+1j.*D.*beta_2.*cos_theta_2), -exp(-1j.*D.*beta_3.*cos_theta_3);
             cos(theta_1)./Z_1, cos_theta_2./Z_2, -cos_theta_2./Z_2, 0;
             0, -(cos_theta_2./Z_2).*exp(-D/2.*alpha_2-1j.*D.*beta_2.*cos_theta_2), (cos_theta_2./Z_2).*exp(D/2.*alpha_2+1j.*D.*beta_2.*cos_theta_2), (cos_theta_3./Z_3).*exp(-1j.*D.*beta_3.*cos_theta_3)];
        TM_A = [1, -1, -1, 0;
             0, exp(-D/2.*alpha_2-1j.*D.*beta_2.*cos_theta_2), exp(D/2.*alpha_2+1j.*D.*beta_2.*cos_theta_2), -exp(-1j.*D.*beta_3.*cos_theta_3);
             cos(theta_1).*Z_1, cos_theta_2.*Z_2, -cos_theta_2.*Z_2, 0;
             0, -(cos_theta_2.*Z_2).*exp(-D/2.*alpha_2-1j.*D.*beta_2.*cos_theta_2), (cos_theta_2.*Z_2).*exp(D/2.*alpha_2+1j.*D.*beta_2.*cos_theta_2), (cos_theta_3.*Z_3).*exp(-1j.*D.*beta_3.*cos_theta_3)];
        TE_B = [-1; 0; cos(theta_1)/Z_1; 0];
        TM_B = [-1; 0; cos(theta_1)*Z_1; 0];
        TE_X = linsolve(TE_A, TE_B);
        TM_X = linsolve(TM_A, TM_B);
        Gamma_TE(j) = TE_X(1);
        Gamma_TM(j) = TM_X(1);
        T_TE(j) = abs(Z_1./Z_3.*abs(TE_X(4))^2);
        T_TM(j) = abs(Z_3./Z_1.*abs(TM_X(4))^2);
    end
    figure(transmitted_TE);
    plot(omega/(2*pi*1e12), 10.*log10(T_TE)); hold on;
    ylim([-2, 0]);
    figure(reflected_TE);
    plot(omega/(2*pi*1e12), 20.*log10(abs(Gamma_TE))); hold on;
    ylim([-10, 0]);
    figure(transmitted_TM);
    plot(omega/(2*pi*1e12), 10.*log10(T_TM)); hold on;
    ylim([-2, 0]);
    figure(reflected_TM);
    plot(omega/(2*pi*1e12), 10.*log10(abs(Gamma_TM))); hold on;
    ylim([-10 0]);
end
figure(transmitted_TE);
title("TE Oblique Incidence: Transmitted power vs. wavelength for 3-layer stack");
xlabel("wavelength [nm]");
ylabel("power delivered to Si (T_{TE}) [dB]");
legend(legend_title);
figure(reflected_TE);
title("TE Oblique Incidence: Reflected power vs. wavelength for 3-layer stack");
xlabel("wavelength [nm]");
ylabel("reflected power (R_{TE}) [dB]");
legend(legend_title);
figure(transmitted_TM);
title("TM Oblique Incidence: Transmitted power vs. wavelength for 3-layer stack");
xlabel("wavelength [nm]");
ylabel("power delivered to Si (T_{TM}) [dB]");
legend(legend_title);
figure(reflected_TM);
title("TM Oblique Incidence: Reflected power vs. wavelength for 3-layer stack");
xlabel("wavelength [nm]");
ylabel("reflected power (R_{TM}) [dB]");
legend(legend_title);