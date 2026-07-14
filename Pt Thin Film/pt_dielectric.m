%analytical dielectric model - pt

clear;

%dielectric function from Drude-Lorentz model
%number of elements
n = 500;

%constant definitions
eps0 = 8.854e-12;   %vacuum permittivity
c = 3e8;   %speed of light
h = 6.626e-34;   %planck's constant (J*s)
hbar = h/(2*pi);   %reduced planck's constant (J*s)
ev_to_j = 1.602e-19;   %eV to J conversion factor
epsinf = 1.76;   %dielectric at infinite frequency

%Drude parameters
h_Omegap_ev = 7.26;   %plasma frequency (eV)
Omegap = h_Omegap_ev * ev_to_j / hbar;

h_Gamma_ev = 0.046;   %relaxation constant (eV)
Gamma = h_Gamma_ev * ev_to_j / hbar;

%Lorentz parameters
omegaj = [0.831, 1.471, 3.449, 10.461].*ev_to_j./hbar;   %resonance frequency
Omegaj = [6.761, 6.895, 9.452, 15.473].*ev_to_j./hbar;   %strength
gammaj = [0.879, 2.052, 5.424, 6.286].*ev_to_j./hbar;   %relaxation constant

%independent variables
thetai = linspace(70, 90, n);   %incident angle (in degrees) on x-axis, linear scale
lambda = linspace(0.19, 3.5, n);   %wavelengths 0.19 to 3.5 μm -> converted to m later
[T,L] = meshgrid(thetai, lambda);   %plotting coordinates

%containers
epsf = zeros([1 n]);   %dielectric function of film

for i = 1:n
    lambdai = lambda(1,i) * 1e-6;   %convert μm to m
    omega = 2*pi*c/lambdai;
    
    %Lorentz oscillators
    lo = 0;
    for j = 1:4
        lo = lo + Omegaj(1,j)^2 / (omegaj(1,j)^2 - omega^2 - 1i*gammaj(1,j)*omega);
    end

    %Drude-Lorentz
    epsf(1,i) = epsinf - Omegap^2 / (omega^2 + 1i*Gamma*omega) + lo;
end


%dielectric function from experimental data
%retrieve experimental data (tselin 80nm)
fileID = fopen('Datasets/Tselin-80nm-n.txt');
tselin_80nm_n = fscanf(fileID, '%f%f', [2 2000]);
fclose(fileID);

fileID = fopen('Datasets/Tselin-80nm-k.txt');
tselin_80nm_k = fscanf(fileID, '%f%f', [2 2000]);
fclose(fileID);

exp_eps = zeros([1 2000]);   %dielectric constant container

%convert refractive indices to dielectric constants
for i = 1:2000
    lambdai = tselin_80nm_n(1,i) * 1e-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;

    N = tselin_80nm_n(2,i) + 1i*tselin_80nm_k(2,i);   %N = n + ik
    eps = N^2;
    exp_eps(1,i) = eps;
end


%compare model to experimental data
clf;
tiledlayout(1,2)

ax1 = nexttile;
hold(ax1, 'on')
plot(ax1, lambda, real(epsf))
plot(ax1, tselin_80nm_n(1,:), real(exp_eps))
hold(ax1, 'off')
title('Real part of dielectric vs. wavelength for 80nm film')
xlabel('Wavelength (μm)')
ylabel('Dielectric function, real part')

ax2 = nexttile;
hold(ax2, 'on')
plot(ax2, lambda, imag(epsf))
plot(ax2, tselin_80nm_n(1,:), imag(exp_eps))
hold(ax2, 'off')
title('Imaginary part of dielectric vs. wavelength for 80nm film')
xlabel('Wavelength (μm)')
ylabel('Dielectric function, imaginary part')
legend('Model','Experimental',Location='northeastoutside')