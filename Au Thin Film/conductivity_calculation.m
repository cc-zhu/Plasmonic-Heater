%computing conductivity from dielectric function (model verification)

clear;

%number of elements
n = 1000;

%constant definitions
eps0 = 8.854e-12;   %vacuum permittivity
c = 3e8;   %speed of light
h = 6.626e-34;   %planck's constant (J*s)
hbar = h/(2*pi);   %reduced planck's constant (J*s)

%independent variables
df = linspace(20, 180, n);   %film thickness (nm), linear scale -> converted to m later
lambda0 = linspace(0.8, 3000, n);   %operating free space wavelength (μm), logarithmic scale -> converted to m later

%fixed parameters
lambda1 = 900e-9;   %fixed wavelength (m) for thickness sweep
omega1 = 2*pi*c/lambda1;   %corresponding omega

%fixed film thicknesses (m)
df25 = 25e-9;
df50 = 50e-9;
df100 = 100e-9;
df150 = 150e-9;
df200 = 200e-9;

%Drude parameters
epsinf = 1.63;   %dielectric constant at infinite frequency

h_omegap_ev = 8.5;   %hbar * omegap = 8.5 eV
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

%(relative) dielectric function of film - Drude model with thickness correction
epsf = @(t,omega) epsinf - omegap^2/(omega^2 + gamma(t)^2) + 1i*gamma(t)*omegap^2/(omega*(omega^2 + gamma(t)^2));   %separated
%epsf = @(t,omega) epsinf - omegap^2/(omega^2 + 1i*gamma(t)*omega);   %compact

%containers - conductivity
s25 = zeros([1 n]);
s50 = zeros([1 n]);
s100 = zeros([1 n]);
s150 = zeros([1 n]);
s200 = zeros([1 n]);

%optical conductivity vs wavelength
for i = 1:n
    %wavelength variation
    lambdai = lambda0(1,i)*10^-6;   %convert μm to m
    omegai = 2*pi*c/(lambdai);

    %25nm film
    dfi = df25;
    eps = epsf(dfi,omegai);   %compute dielectric function
    s25(1,i) = -1i*omegai*eps0*(eps - 1);   %compute (complex) conductivity
    
    %50nm film
    dfi = df50;
    eps = epsf(dfi,omegai);
    s50(1,i) = -1i*omegai*eps0*(eps - 1);

    %100nm film
    dfi = df100;
    eps = epsf(dfi,omegai);
    s100(1,i) = -1i*omegai*eps0*(eps - 1);
    
    %150nm film
    dfi = df150;
    eps = epsf(dfi,omegai);
    s150(1,i) = -1i*omegai*eps0*(eps - 1);

    %200nm film
    dfi = df200;
    eps = epsf(dfi,omegai);
    s200(1,i) = -1i*omegai*eps0*(eps - 1);
end

clf
tiledlayout(2,1)

ax1 = nexttile;
hold(ax1,'on')
plot(ax1, lambda0, real(s25))
plot(ax1, lambda0, real(s50))
plot(ax1, lambda0, real(s100))
plot(ax1, lambda0, real(s150))
plot(ax1, lambda0, real(s200))
hold(ax1,'off')
title('Real part of conductivity vs. wavelength')
xlabel('Wavelength (μm)')
ylabel('Conductivity (S/m)')
lgd = legend('25 nm','50 nm','100 nm','150 nm','200 nm',Location='northeastoutside');
title(lgd,'Film thicknesses')

ax2 = nexttile;
hold(ax2,'on')
plot(ax2, lambda0, imag(s25))
plot(ax2, lambda0, imag(s50))
plot(ax2, lambda0, imag(s100))
plot(ax2, lambda0, imag(s150))
plot(ax2, lambda0, imag(s200))
hold(ax2,'off')
title('Imaginary part of conductivity vs. wavelength')
xlabel('Wavelength (μm)')
ylabel('Conductivity (S/m)')