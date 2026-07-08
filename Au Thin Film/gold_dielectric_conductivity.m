%exploring behavior of conductivity and dielectric function from refractive index data

clear;

%retrieve experimental data (yakubovsky 25nm, 53nm, 117nm)
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
df25 = 25e-9;
df53 = 53e-9;
df117 = 117e-9;

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

%(relative) dielectric function of film - Drude model with thickness correction
epsf = @(t,omega) epsinf - omegap^2/(omega^2 + gamma(t)^2) + 1i*gamma(t)*omegap^2/(omega*(omega^2 + gamma(t)^2));


%containers
%dielectric constant from experimental refractive index
exp_eps_25nm = zeros([1 171]);
exp_eps_53nm = zeros([1 171]);
exp_eps_117nm = zeros([1 171]);

%conductivity
exp_s_25nm = zeros([1 171]);
exp_s_53nm = zeros([1 171]);
exp_s_117nm = zeros([1 171]);

%dielectric constants from analytical model
model_eps_25nm = zeros([1 n]);
model_eps_53nm = zeros([1 n]);
model_eps_117nm = zeros([1 n]);

%conductivity
model_s_25nm = zeros([1 n]);
model_s_53nm = zeros([1 n]);
model_s_117nm = zeros([1 n]);


%convert refractive indices to dielectric constants & calculate conductivity
for i = 1:171
    %25 nm film
    lambdai = yakubovsky_25nm_n(1,i)*10^-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;

    N = yakubovsky_25nm_n(2,i) + 1i*yakubovsky_25nm_k(2,i);   %N = n + ik
    eps = N^2;
    exp_eps_25nm(1,i) = eps;
    exp_s_25nm(1,i) = -1i*omegai*eps0*(eps - 1);

    %53 nm film
    lambdai = yakubovsky_53nm_n(1,i)*10^-6;
    omegai = 2*pi*c/lambdai;

    N = yakubovsky_53nm_n(2,i) + 1i*yakubovsky_53nm_k(2,i);
    eps = N^2;
    exp_eps_53nm(1,i) = eps;
    exp_s_53nm(1,i) = -1i*omegai*eps0*(eps - 1);

    %117 nm film
    lambdai = yakubovsky_117nm_n(1,i)*10^-6;
    omegai = 2*pi*c/lambdai;

    N = yakubovsky_117nm_n(2,i) + 1i*yakubovsky_117nm_k(2,i);
    eps = N^2;
    exp_eps_117nm(1,i) = eps;
    exp_s_117nm(1,i) = -1i*omegai*eps0*(eps - 1);
end


%calculate dielectric constants & conductivity from model
for i = 1:n
    lambdai = lambda0(1,i)*10^-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;
    
    %25nm film
    eps = epsf(df25,omegai);
    model_eps_25nm(1,i) = eps;
    model_s_25nm(1,i) = -1i*omegai*eps0*(eps - 1);

    %53nm film
    eps = epsf(df53,omegai);
    model_eps_53nm(1,i) = eps;
    model_s_53nm(1,i) = -1i*omegai*eps0*(eps - 1);

    %117nm film
    eps = epsf(df117,omegai);
    model_eps_117nm(1,i) = eps;
    model_s_117nm(1,i) = -1i*omegai*eps0*(eps - 1);
end

% %plots - dielectric function from experimental data & analytical model
% clf
% tiledlayout(3,2)
% 
% %25nm film
% ax1 = nexttile;
% hold(ax1,'on')
% plot(ax1, yakubovsky_25nm_n(1,:), real(exp_eps_25nm))
% plot(ax1, lambda0, real(model_eps_25nm))
% hold(ax1,'off')
% title('Real part of dielectric vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectric function, real part')
% 
% ax2 = nexttile;
% hold(ax2,'on')
% plot(ax2, yakubovsky_25nm_n(1,:), imag(exp_eps_25nm))
% plot(ax2, lambda0, imag(model_eps_25nm))
% hold(ax2,'off')
% title('Imaginary part of dielectric vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectric function, imaginary part')
% legend('Experimental', 'Model', Location='northeastoutside')
% 
% 
% %53nm film
% ax3 = nexttile;
% hold(ax3,'on')
% plot(ax3, yakubovsky_53nm_n(1,:), real(exp_eps_53nm))
% plot(ax3, lambda0, real(model_eps_53nm))
% hold(ax3,'off')
% title('Real part of dielectric vs. wavelength for 53nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectric function, real part')
% 
% ax4 = nexttile;
% hold(ax4,'on')
% plot(ax4, yakubovsky_53nm_n(1,:), imag(exp_eps_53nm))
% plot(ax4, lambda0, imag(model_eps_53nm))
% hold(ax4,'off')
% title('Imaginary part of dielectric vs. wavelength for 53nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectric function, imaginary part')
% 
% 
% %117nm film
% ax5 = nexttile;
% hold(ax5,'on')
% plot(ax5, yakubovsky_117nm_n(1,:), real(exp_eps_117nm))
% plot(ax5, lambda0, real(model_eps_117nm))
% hold(ax5,'off')
% title('Real part of dielectric vs. wavelength for 117nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectric function, real part')
% 
% ax6 = nexttile;
% hold(ax6,'on')
% plot(ax6, yakubovsky_117nm_n(1,:), imag(exp_eps_117nm))
% plot(ax6, lambda0, imag(model_eps_117nm))
% hold(ax6,'off')
% title('Imaginary part of dielectric vs. wavelength for 117nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectic function, imaginary part')

% %conductivity vs. dielectric relation
% clf
% tiledlayout(2,1)
% 
% ax1 = nexttile;
% hold(ax1,'on')
% plot(ax1, imag(model_eps_25nm), real(model_s_25nm))
% plot(ax1, imag(model_eps_53nm), real(model_s_53nm))
% plot(ax1, imag(model_eps_117nm), real(model_s_117nm))
% hold(ax1,'off')
% title('Real conductivity vs. imaginary dielectric')
% xlabel('Imaginary part of dielectric function')
% ylabel('Real part of conductivity (S/m)')
% xlim(ax1, [0 20])
% lgd = legend('25 nm', '53 nm', '117 nm',Location='northeastoutside');
% title(lgd,'Film thicknesses')
% 
% ax2 = nexttile;
% hold(ax2,'on')
% plot(ax2, real(model_eps_25nm), imag(model_s_25nm))
% plot(ax2, real(model_eps_53nm), imag(model_s_53nm))
% plot(ax2, real(model_eps_117nm), imag(model_s_117nm))
% hold(ax2,'off')
% title('Imaginary conductivity vs. real dielectric')
% xlabel('Real part of dielectric function')
% ylabel('Imaginary part of conductivity (S/m)')
% xlim(ax2, [-200 0])


%dielectric required for 10^6 S/m conductivity
lambda = linspace(0.3, 10, n);   %wavelengths 0.3 to 10 μm -> converted to m later
sigma = 1e6;   %target conductivity
req_eps_i = zeros([1 n]);   %required dielectric function, imaginary part
model_eps_i = zeros([1 n]);   %dielectric function, imaginary part given by model
df1 = 100e-9;   %100 nm film

for i = 1:n
    lambdai = lambda(1,i)*10^-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;
    
    req_eps = 1 + 1i*sigma/(omegai*eps0);
    req_eps_i(1,i) = imag(req_eps);

    eps = epsf(df1,omegai);
    model_eps_i(1,i) = imag(eps);
end

%plot required imaginary dielectric and imaginary dielectric given by model
clf
hold on
plot(lambda, req_eps_i)
plot(lambda, model_eps_i)
hold off
title('Imaginary dielectric vs. wavelength for 100 nm film')
xlabel('Wavelength (μm)')
ylabel('Dielectric function, imaginary part')
legend('Required for 1e6 S/m conductivity','Model',Location='northeastoutside')