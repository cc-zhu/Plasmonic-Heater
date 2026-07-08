%exploring behavior of conductivity and dielectric function from refractive index data

clear;

%retrieve experimental data (tselin 80nm)
fileID = fopen('Datasets/Tselin-80nm-n.txt');
tselin_80nm_n = fscanf(fileID, '%f%f', [2 2000]);
fclose(fileID);

fileID = fopen('Datasets/Tselin-80nm-k.txt');
tselin_80nm_k = fscanf(fileID, '%f%f', [2 2000]);
fclose(fileID);

%containers - values calculated from experimental refractive index
exp_eps = zeros([1 2000]);   %dielectric constant
exp_s = zeros([1 2000]);   %conductivity

%constants
eps0 = 8.854e-12;   %vacuum permittivity
c = 3e8;   %speed of light


%convert refractive indices to dielectric constants & calculate conductivity
for i = 1:2000
    lambdai = tselin_80nm_n(1,i) * 1e-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;

    N = tselin_80nm_n(2,i) + 1i*tselin_80nm_k(2,i);   %N = n + ik
    eps = N^2;
    exp_eps(1,i) = eps;
    exp_s(1,i) = -1i*omegai*eps0*(eps - 1);
end


%plot conductivity vs. wavelength
clf
hold on
plot(tselin_80nm_n(1,:), real(exp_s))
plot(tselin_80nm_n(1,:), imag(exp_s))
hold off
title('Conductivity vs. wavelength for 80nm film')
xlabel('Wavelength (μm)')
ylabel('Conductivity (S/m)')
legend('Real part', 'Imaginary part', Location='northeastoutside')


%number of elements
n = 1000;

%dielectric required for 10^6 S/m conductivity
lambda = linspace(0.19, 3.5, n);   %wavelengths 0.19 to 3.5 μm -> converted to m later
sigma = 1e6;   %target conductivity
req_eps_i = zeros([1 n]);   %required dielectric function, imaginary part
df1 = 80e-9;   %100 nm film

for i = 1:n
    lambdai = lambda(1,i) * 1e-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;
    
    req_eps = 1 + 1i*sigma/(omegai*eps0);
    req_eps_i(1,i) = imag(req_eps);
end

% %plot required imaginary dielectric and imaginary dielectric given by model
% clf
% hold on
% plot(lambda, req_eps_i)
% plot(tselin_80nm_n(1,:), imag(exp_eps))
% hold off
% title('Imaginary dielectric vs. wavelength for 80 nm film')
% xlabel('Wavelength (μm)')
% ylabel('Dielectric function, imaginary part')
% legend('Required for 1e6 S/m conductivity','Experimental',Location='northeastoutside')