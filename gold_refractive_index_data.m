%verifying analytical model by plotting refractive indices and conductivities against experimental data

clear;

%retrieve experimental data (rosenblatt 11nm, 21nm, 44nm; yakubovsky 25nm, 53nm, 117nm)

%rosenblatt datasets
% 11 nm
fileID = fopen('Datasets/Rosenblatt-11nm-n.txt');
rosenblatt_11nm_n = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

fileID = fopen('Datasets/Rosenblatt-11nm-k.txt');
rosenblatt_11nm_k = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

% 21 nm
fileID = fopen('Datasets/Rosenblatt-21nm-n.txt');
rosenblatt_21nm_n = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

fildID = fopen('Datasets/Rosenblatt-21nm-k.txt');
rosenblatt_21nm_k = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

% 44 nm
fileID = fopen('Datasets/Rosenblatt-44nm-n.txt');
rosenblatt_44nm_n = fscanf(fileID, '%f%f', [2 176]);
fclose(fileID);

fileID = fopen('Datasets/Rosenblatt-44nm-k.txt');
rosenblatt_44nm_k = fscanf(fileID, '%f%f', [2 176]);
fclose(fileID);

%yakubovsky datasets
% 25 nm
fileID = fopen('Datasets/Yakubovsky-25nm-n.txt');
yakubovsky_25nm_n = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

fileID = fopen('Datasets/Yakubovsky-25nm-k.txt');
yakubovsky_25nm_k = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

% 53 nm
fileID = fopen('Datasets/Yakubovsky-53nm-n.txt');
yakubovsky_53nm_n = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

fileID = fopen('Datasets/Yakubovsky-53nm-k.txt');
yakubovsky_53nm_k = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

% 117 nm
fileID = fopen('Datasets/Yakubovsky-117nm-n.txt');
yakubovsky_117nm_n = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

fileID = fopen('Datasets/Yakubovsky-117nm-k.txt');
yakubovsky_117nm_k = fscanf(fileID, '%f%f', [2 171]);
fclose(fileID);

%analytical model
%number of elements
n = 500;

%constant definitions
eps0 = 8.854e-12;   %vacuum permittivity
c = 3e8;   %speed of light
h = 6.626e-34;   %planck's constant (J*s)
hbar = h/(2*pi);   %reduced planck's constant (J*s)

%independent variables
lambda0 = linspace(0.3, 2, n);   %operating free space wavelength (μm), linear scale -> converted to m later

%fixed thicknesses (m)
df9 = 9e-9;
df25 = 25e-9;
df53 = 53e-9;
df117 = 117e-9;

df11 = 11e-9;
df21 = 21e-9;
df44 = 44e-9;

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
%F = @(t) (1/(1 - 1.5*alpha(t) + 3*alpha(t)^2 - 3*alpha(t)^3*log(1 + 1/alpha(t))));

gamma = @(t) gammaep + gammagb(t) + gammas;   %total damping parameter (1/s)
%gamma = @(t) gammaep*F(t) + gammas;
%gamma = @(t) gammaep;   %without thickness correction

%(relative) dielectric function of film - Drude model with thickness correction
epsf = @(t,omega) epsinf - omegap^2/(omega^2 + gamma(t)^2) + 1i*gamma(t)*omegap^2/(omega*(omega^2 + gamma(t)^2));   %separated
%epsf = @(t,omega) epsinf - omegap^2/(omega^2 + 1i*gamma(t)*omega);   %compact

%containers
%refractive index, real part
n25 = zeros([1 n]);   %25nm
n53 = zeros([1 n]);   %53nm
n117 = zeros([1 n]);   %117nm

n11 = zeros([1 n]);   %11nm
n21 = zeros([1 n]);   %21nm
n44 = zeros([1 n]);   %44nm

%refractive index, imaginary part
k25 = zeros([1 n]);
k53 = zeros([1 n]);
k117 = zeros([1 n]);

k11 = zeros([1 n]);
k21 = zeros([1 n]);
k44 = zeros([1 n]);

%conductivity, real part
sr25 = zeros([1 n]);
sr53 = zeros([1 n]);
sr117 = zeros([1 n]);

%conductivity, imaginary part
si25 = zeros([1 n]);
si53 = zeros([1 n]);
si117 = zeros([1 n]);

for i = 1:n
    lambdai = lambda0(1,i)*10^-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;
    
    %yakubovsky
    %25nm film
    eps = epsf(df25,omegai);
    n = sqrt(eps);
    n25(1,i) = real(n);
    k25(1,i) = imag(n);

    sigma = -1i*omegai*eps0*(eps - 1);
    sr25(1,i) = real(sigma);
    si25(1,i) = imag(sigma);

    %53nm film
    eps = epsf(df53,omegai);
    n = sqrt(eps);
    n53(1,i) = real(n);
    k53(1,i) = imag(n);

    sigma = -1i*omegai*eps0*(eps - 1);
    sr53(1,i) = real(sigma);
    si53(1,i) = imag(sigma);

    %117nm film
    eps = epsf(df117,omegai);
    n = sqrt(eps);
    n117(1,i) = real(n);
    k117(1,i) = imag(n);

    sigma = -1i*omegai*eps0*(eps - 1);
    sr117(1,i) = real(sigma);
    si117(1,i) = imag(sigma);


    %rosenblatt
    %11nm film
    eps = epsf(df11,omegai);
    n = sqrt(eps);
    n11(1,i) = real(n);
    k11(1,i) = imag(n);

    %21nm film
    eps = epsf(df21,omegai);
    n = sqrt(eps);
    n21(1,i) = real(n);
    k21(1,i) = imag(n);

    %44nm film
    eps = epsf(df44,omegai);
    n = sqrt(eps);
    n44(1,i) = real(n);
    k44(1,i) = imag(n);
end

% %plots - refractive indices from experimental data & analytical model
% clf
% tiledlayout(3,2)   %yakubovsky
% 
% %yakubovsky 25nm
% ax1 = nexttile;
% hold(ax1,'on')
% plot(ax1, yakubovsky_25nm_n(1,:), yakubovsky_25nm_n(2,:))
% plot(ax1, lambda0, n25)
% hold(ax1,'off')
% title('Real part of refractive index vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, real part')
% 
% ax2 = nexttile;
% hold(ax2,'on')
% plot(ax2, yakubovsky_25nm_k(1,:), yakubovsky_25nm_k(2,:))
% plot(ax2, lambda0, k25)
% hold(ax2,'off')
% title('Imaginary part of refractive index vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, imaginary part')
% legend('Experimental data', 'Model', Location='northeastoutside')
% 
% 
% %yakubovsky 53nm
% ax3 = nexttile;
% hold(ax3,'on')
% plot(ax3, yakubovsky_53nm_n(1,:), yakubovsky_53nm_n(2,:))
% plot(ax3, lambda0, n53)
% hold(ax3,'off')
% title('Real part of refractive index vs. wavelength for 53nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, real part')
% 
% ax4 = nexttile;
% hold(ax4,'on')
% plot(ax4, yakubovsky_53nm_k(1,:), yakubovsky_53nm_k(2,:))
% plot(ax4, lambda0, k53)
% hold(ax4,'off')
% title('Imaginary part of refractive index vs. wavelength for 53nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, imaginary part')
% 
% 
% %yakubovsky 117nm
% ax5 = nexttile;
% hold(ax5,'on')
% plot(ax5, yakubovsky_117nm_n(1,:), yakubovsky_117nm_n(2,:))
% plot(ax5, lambda0, n117)
% hold(ax5,'off')
% title('Real part of refractive index vs. wavelength for 117nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, real part')
% 
% ax6 = nexttile;
% hold(ax6,'on')
% plot(ax6, yakubovsky_117nm_k(1,:), yakubovsky_117nm_k(2,:))
% plot(ax6, lambda0, k117)
% hold(ax6,'off')
% title('Imaginary part of refractive index vs. wavelength for 117nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, imaginary part')


% tiledlayout(3,2)   %rosenblatt
% 
% %rosenblatt 11nm
% ax1 = nexttile;
% hold(ax1,'on')
% plot(ax1, rosenblatt_11nm_n(1,:), rosenblatt_11nm_n(2,:))
% plot(ax1, lambda0, n11)
% hold(ax1,'off')
% title('Real part of refractive index vs. wavelength for 11nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, real part')
% 
% ax2 = nexttile;
% hold(ax2,'on')
% plot(ax2, rosenblatt_11nm_k(1,:), rosenblatt_11nm_k(2,:))
% plot(ax2, lambda0, k11)
% hold(ax2,'off')
% title('Imaginary part of refractive index vs. wavelength for 11nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, imaginary part')
% legend('Experimental data', 'Model', Location='northeastoutside')
% 
% 
% %rosenblatt 21nm
% ax3 = nexttile;
% hold(ax3,'on')
% plot(ax3, rosenblatt_21nm_n(1,:), rosenblatt_21nm_n(2,:))
% plot(ax3, lambda0, n21)
% hold(ax3,'off')
% title('Real part of refractive index vs. wavelength for 21nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, real part')
% 
% ax4 = nexttile;
% hold(ax4,'on')
% plot(ax4, rosenblatt_21nm_k(1,:), rosenblatt_21nm_k(2,:))
% plot(ax4, lambda0, k21)
% hold(ax4,'off')
% title('Imaginary part of refractive index vs. wavelength for 21nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, imaginary part')
% 
% 
% %rosenblatt 44nm
% ax5 = nexttile;
% hold(ax5,'on')
% plot(ax5, rosenblatt_44nm_n(1,:), rosenblatt_44nm_n(2,:))
% plot(ax5, lambda0, n44)
% hold(ax5,'off')
% title('Real part of refractive index vs. wavelength for 44nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, real part')
% 
% ax6 = nexttile;
% hold(ax6,'on')
% plot(ax6, rosenblatt_44nm_k(1,:), rosenblatt_44nm_k(2,:))
% plot(ax6, lambda0, k44)
% hold(ax6,'off')
% title('Imaginary part of refractive index vs. wavelength for 44nm film')
% xlabel('Wavelength (μm)')
% ylabel('Refractive index, imaginary part')


%calculate conductivity from dielectric data
s25 = zeros([1 171]);
s53 = zeros([1 171]);
s117 = zeros([1 171]);

for i = 1:171
    %yakubovsky
    %25nm film
    ni = yakubovsky_25nm_n(2,i);
    ki = yakubovsky_25nm_k(2,i);
    lambdai = yakubovsky_25nm_n(1,i)*10^-6;
    omegai = 2*pi*c/lambdai;

    epsi = (ni + 1i*ki)^2;
    s25(1,i) = -1i*omegai*eps0*(epsi - 1);

    %53nm film
    ni = yakubovsky_53nm_n(2,i);
    ki = yakubovsky_53nm_k(2,i);
    lambdai = yakubovsky_53nm_n(1,i)*10^-6;
    omegai = 2*pi*c/lambdai;

    epsi = (ni + 1i*ki)^2;
    s53(1,i) = -1i*omegai*eps0*(epsi - 1);

    %117nm film
    ni = yakubovsky_117nm_n(2,i);
    ki = yakubovsky_117nm_k(2,i);
    lambdai = yakubovsky_117nm_n(1,i)*10^-6;
    omegai = 2*pi*c/lambdai;

    epsi = (ni + 1i*ki)^2;
    s117(1,i) = -1i*omegai*eps0*(epsi - 1);
end

% clf
% tiledlayout(3,2)   %yakubovsky
% 
% %yakubovsky 25nm
% ax1 = nexttile;
% hold(ax1,'on')
% plot(ax1, yakubovsky_25nm_n(1,:), real(s25))
% plot(ax1, lambda0, sr25)
% hold(ax1,'off')
% title('Real part of conductivity vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')
% 
% ax2 = nexttile;
% hold(ax2,'on')
% plot(ax2, yakubovsky_25nm_n(1,:), imag(s25))
% plot(ax2, lambda0, si25)
% hold(ax2,'off')
% title('Imaginary part of conductivity vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')
% legend('Experimental data', 'Model', Location='northeastoutside')
% 
% %yakubovsky 53nm
% ax3 = nexttile;
% hold(ax3,'on')
% plot(ax3, yakubovsky_53nm_n(1,:), real(s53))
% plot(ax3, lambda0, sr53)
% hold(ax3,'off')
% title('Real part of conductivity vs. wavelength for 53nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')
% 
% ax4 = nexttile;
% hold(ax4,'on')
% plot(ax4, yakubovsky_53nm_n(1,:), imag(s53))
% plot(ax4, lambda0, si53)
% hold(ax4,'off')
% title('Imaginary part of conductivity vs. wavelength for 53nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')
% 
% %yakubovsky 117nm
% ax5 = nexttile;
% hold(ax5,'on')
% plot(ax5, yakubovsky_117nm_n(1,:), real(s117))
% plot(ax5, lambda0, sr117)
% hold(ax5,'off')
% title('Real part of conductivity vs. wavelength for 117nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')
% 
% ax6 = nexttile;
% hold(ax6,'on')
% plot(ax6, yakubovsky_117nm_n(1,:), imag(s117))
% plot(ax6, lambda0, si117)
% hold(ax6,'off')
% title('Imaginary part of conductivity vs. wavelength for 117nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')

