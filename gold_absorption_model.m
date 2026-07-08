%analytical absorption model - gold

clear;

%number of elements
n = 500;

%constant definitions
eps0 = 8.854e-12;   %vacuum permittivity
c = 3e8;   %speed of light
h = 6.626e-34;   %planck's constant (J*s)
hbar = h/(2*pi);   %reduced planck's constant (J*s)

%independent variables
thetai = linspace(70, 90, n);   %incident angle (degrees), linear scale -> converted to rad later
df = linspace(20, 200, n);   %film thickness (nm), linear scale -> converted to m later
lambda0 = linspace(100, 8e4, n);   %operating free space wavelength (nm), linear scale -> converted to m later

%fixed parameters
lambda1 = 60e-6;   %fixed wavelength (m) for thickness sweep
omega1 = 2*pi*c/lambda1;   %corresponding omega
df1 = 100e-9;   %fixed thickness (m) for wavelength sweep

[T1, Df] = meshgrid(thetai, df);   %plotting coordinates for varying thickness
[T2, L0] = meshgrid(thetai, lambda0);   %plotting coordinates for varying wavelength
A1 = zeros(n); A2 = zeros(n);   %containers for absorption values

%Drude parameters
epsinf = 1.63;   %dielectric constant at infinite frequency

h_omegap_ev = 8.87;   %hbar * omegap = 8.5 eV
h_omegap = h_omegap_ev * 1.602e-19;   %convert hbar * omegap from eV to J
omegap = h_omegap/hbar;   %plasma frequency (1/s)

%damping rate
gammaep = 4.6e13;   %electron-phonon scattering (1/s)
gammas = 1.35e13;   %scattering on film surface (1/s)

R = 0.31;   %grain-boundary reflection coefficient
Lambda = 39e-9;   %mean free path for conducting electrons in bulk gold (m)
D = @(t) (1.522e-8)*log(t) + 2.847e-7;   %average crystalite size as a function of thickness (m)
alpha = @(t) Lambda*R/(D(t)*(1 - R));

%electron-grain boundary scattering rate (1/s)
gammagb = @(t) gammaep*(1/(1 - 1.5*alpha(t) + 3*alpha(t)^2 - 3*alpha(t)^3*log(1 + 1/alpha(t))) - 1);

gamma = @(t) gammaep + gammagb(t) + gammas;   %total damping parameter (1/s)

%dielectric function of film - Drude model with thickness correction
epsf = @(t,omega) epsinf - omegap^2/(omega^2 + gamma(t)^2) + 1i*gamma(t)*omegap^2/(omega*(omega^2 + gamma(t)^2));

%transfer matrix method
k0 = @(omega) omega/c;   %free space wave number (1/m)

% %simplified transfer matrix terms
% M11 = @(t,omega) cos(sqrt(epsf(t,omega)) * k0(omega) * t);
% M12 = @(t,omega) 1i*sin(sqrt(epsf(t,omega))* k0(omega) * t)/sqrt(epsf(t,omega));
% M21 = @(t,omega) 1i*sqrt(epsf(t,omega))*sin(sqrt(epsf(t,omega)) * k0(omega) * t);
% M22 = @(t,omega) cos(sqrt(epsf(t,omega)) * k0(omega) * t);

%full transfer matrix terms
M11 = @(t,omega,theta) cos(sqrt(epsf(t,omega)-(sin(theta))^2) * k0(omega) * t);
M12 = @(t,omega,theta) 1i*sin(sqrt(epsf(t,omega)-(sin(theta))^2)* k0(omega) * t)*sqrt((epsf(t,omega)-(sin(theta))^2)/(epsf(t,omega))^2);
M21 = @(t,omega,theta) 1i*sin(sqrt(epsf(t,omega)-(sin(theta))^2) * k0(omega) * t)*sqrt(((epsf(t,omega))^2)/(epsf(t,omega)-(sin(theta))^2));
M22 = @(t,omega,theta) cos(sqrt(epsf(t,omega)-(sin(theta))^2) * k0(omega) * t);

eta = @(theta) 1/cos(theta);

%reflection coefficient
r = @(t,omega,theta) (eta(theta)*M22(t,omega,theta) - eta(theta)*M11(t,omega,theta) - ((eta(theta)^2)*M12(t,omega,theta) - M21(t,omega,theta)))/(eta(theta)*M22(t,omega,theta) + eta(theta)*M11(t,omega,theta) - ((eta(theta)^2)*M12(t,omega,theta) + M21(t,omega,theta)));
%transmission coefficient
t = @(t,omega,theta) 2*eta(theta)/(eta(theta)*M22(t,omega,theta) + eta(theta)*M11(t,omega,theta) - ((eta(theta)^2)*M12(t,omega,theta) + M21(t,omega,theta)));


% %check general model behavior
% eps = zeros([1 n]);
% eps1 = zeros([1 n]);
% eps2 = zeros([1 n]);
% g = zeros([1 n]);
% Dt = zeros([1 n]);
% 
% for i = 1:n
%     % %thickness variation
%     % dfi = df(1,i)*10^-9;   %convert nm to m
%     % lambdai = lambda1;
%     % omegai = 2*pi*c/lambdai;
% 
%     %wavelength variation
%     dfi = df1;
%     lambdai = lambda0(1,i)*10^-9;   %convert nm to m
%     omegai = 2*pi*c/lambdai;
% 
%     eps(1,i) = epsf(dfi,omegai);
%     eps1(1,i) = real(eps(1,i));
%     eps2(1,i) = imag(eps(1,i));
% end
% 
% for i = 1:n
%     %thickness variation
%     dfi = df(1,i)*10^-9;   %convert nm to m
%     lambdai = lambda1;
%     omegai = 2*pi*c/lambdai;
% 
%     % %wavelength variation
%     % dfi = df1;
%     % lambdai = lambda0(1,i)*10^-9;   %convert nm to m
%     % omegai = 2*pi*c/lambdai;
% 
%     Di = D(dfi);
%     Dt(1,i) = Di;
% 
%     gi = gamma(dfi);
%     g(1,i) = gi;
% end
% 
% clf
% tiledlayout(2,2)
% 
% ax1 = nexttile;
% plot(ax1,lambda0,eps1)   %real part of dielectric vs. wavelength
% title('Real part of dielectric vs. operating wavelength with df = 40nm')
% xlabel('Operating wavelength (nm)')
% ylabel('Dielectric function, real part')
% 
% ax2 = nexttile;
% plot(ax2,lambda0,eps2)   %imag part of dielectric vs. wavelength
% title('Imaginary part of dielectric vs. operating wavelength with df = 40nm')
% xlabel('Operating wavelength (nm)')
% ylabel('Dielectric function, imaginary part')
% 
% ax3 = nexttile;
% plot(ax3,df,g)   %total damping vs. film thickness
% title('Total damping vs. film thickness with lambda = 900nm')
% xlabel('Film thickness (nm)')
% ylabel('Total damping (1/s)')
% 
% ax4 = nexttile;
% plot(ax4,df,Dt)   %average crystallite size vs. film thickness
% title('Average crystallite size vs. film thickness with lambda = 900nm')
% xlabel('Film thickness (nm)')
% ylabel('Average crystallite size (m)')


%thickness sweep, fixed wavelength
for i = 1:n
    for j = 1:n
        thetaiij = deg2rad(T1(i,j));   %current thetai, convert deg to rad
        dfij = Df(i,j)*10^-9;   %current df, convert nm to m

        rij = r(dfij,omega1,thetaiij);   %reflection coefficient
        tij = t(dfij,omega1,thetaiij);   %transmission coefficient
        Aij = 1 - abs(rij)^2 - abs(tij)^2;   %absorption
        A1(i,j) = Aij;   %assign absorption
    end
end

%max absorption from thickness sweep
Amax_t = max(A1,[],"all")


%wavelength sweep, fixed thickness
for i = 1:n
    for j = 1:n
        thetaiij = deg2rad(T2(i,j));   %current thetai, convert deg to rad
        lambda0ij = L0(i,j)*10^-9;   %current lambda0, convert nm to m
        omegaij = 2*pi*c/lambda0ij;   %corresponding omega

        rij = r(df1,omegaij,thetaiij);   %reflection coefficient
        tij = t(df1,omegaij,thetaiij);   %transmission coefficient
        Aij = 1 - abs(rij)^2 - abs(tij)^2;   %absorption
        A2(i,j) = Aij;   %assign absorption
    end
end

%max absorption from wavelength sweep
Amax_w = max(A2,[],"all")

clf
tiledlayout(1,2)

%plot absorption vs. film thickness and incident angle
ax1 = nexttile;
surf(ax1,T2,Df,A1,'EdgeColor','None');   %3d plot of A1
view(2);   %2d view
colormap hot
colorbar
xlim([80,90])
ylim([20 200])
title('Absorption vs. film thickness and incident angle with lambda = 60 μm')
xlabel('Incident angle (degree)')
ylabel('Film thickness (nm)')

%plot absorption vs. wavelength and incident angle
ax2 = nexttile;
surf(ax2,T1,L0,A2,'EdgeColor','None');   %3d plot of A2
view(2);   %2d view
colormap hot
colorbar
xlim([80,90])
ylim([100 8e4])
title('Absorption vs. operating wavelength and incident angle with df = 100 nm')
xlabel('Incident angle (degree)')
ylabel('Operating wavelength (nm)')