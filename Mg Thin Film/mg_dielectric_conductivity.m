%exploring behavior of conductivity and dielectric function from refractive index data

clear;

%retrieve experimental data (palm 25nm)
fileID = fopen('Datasets/Palm-25nm-n.txt');
palm_25nm_n = fscanf(fileID, '%f%f', [2 670]);
fclose(fileID);

fileID = fopen('Datasets/Palm-25nm-k.txt');
palm_25nm_k = fscanf(fileID, '%f%f', [2 670]);
fclose(fileID);

%containers - values calculated from experimental refractive index
exp_eps = zeros([1 670]);   %dielectric constant
exp_s = zeros([1 670]);   %conductivity

%constants
eps0 = 8.854e-12;   %vacuum permittivity
c = 3e8;   %speed of light


%convert refractive indices to dielectric constants & calculate conductivity
for i = 1:670
    lambdai = palm_25nm_n(1,i) * 1e-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;

    N = palm_25nm_n(2,i) + 1i*palm_25nm_k(2,i);   %N = n + ik
    eps = N^2;
    exp_eps(1,i) = eps;
    exp_s(1,i) = -1i*omegai*eps0*(eps - 1);
end


% %plot conductivity vs. wavelength
% clf
% hold on
% plot(palm_25nm_n(1,:), real(exp_s))
% plot(palm_25nm_n(1,:), imag(exp_s))
% hold off
% title('Conductivity vs. wavelength for 25nm film')
% xlabel('Wavelength (μm)')
% ylabel('Conductivity (S/m)')
% legend('Real part', 'Imaginary part', Location='northeastoutside')


%number of elements
n = 1000;

%dielectric required for 10^6 S/m conductivity
lambda = linspace(0.250019531, 1.684528687, n);   %wavelengths 0.19 to 3.5 μm -> converted to m later
sigma = 1e6;   %target conductivity
req_eps_i = zeros([1 n]);   %required dielectric function, imaginary part
df1 = 25e-9;   %25 nm film

for i = 1:n
    lambdai = lambda(1,i) * 1e-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;
    
    req_eps = 1 + 1i*sigma/(omegai*eps0);
    req_eps_i(1,i) = imag(req_eps);
end

%plot required imaginary dielectric and imaginary dielectric given by model
clf
hold on
plot(lambda, req_eps_i)
plot(palm_25nm_n(1,:), imag(exp_eps))
hold off
title('Imaginary dielectric vs. wavelength for 25 nm film')
xlabel('Wavelength (μm)')
ylabel('Dielectric function, imaginary part')
legend('Required for 1e6 S/m conductivity','Experimental',Location='northeastoutside')