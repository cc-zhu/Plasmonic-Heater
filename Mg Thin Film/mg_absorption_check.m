%absorption of 25nm mg film using experimental refractive index

clear;

%dielectric function from experimental data

%retrieve experimental data (palm 25nm)
fileID = fopen('Datasets/Palm-25nm-n.txt');
palm_25nm_n = fscanf(fileID, '%f%f', [2 670]);
fclose(fileID);

fileID = fopen('Datasets/Palm-25nm-k.txt');
palm_25nm_k = fscanf(fileID, '%f%f', [2 670]);
fclose(fileID);

exp_eps = zeros([1 670]);   %dielectric constant container

%constants
c = 3e8;   %speed of light
ep0 = 8.85419e-12;   %vacuum permittivity
mu0 = 4*pi*1e-7;   %vacuum permeability

%convert refractive indices to dielectric constants
for i = 1:670
    lambdai = palm_25nm_n(1,i) * 1e-6;   %convert μm to m
    omegai = 2*pi*c/lambdai;

    N = palm_25nm_n(2,i) + 1i*palm_25nm_k(2,i);   %N = n + ik
    eps = N^2;
    exp_eps(1,i) = eps;
end


%absorption calculation

%number of elements
n = 670;

%incident angle lower bound
t1 = 55;

%independent variables
thetai = linspace(t1, 90, n);   %incident angle (in degrees) on x-axis, linear scale
lambda = linspace(0.250019531, 1.684528687, n);   %wavelengths 0.19 to 3.5 μm -> converted to m later
[T,L] = meshgrid(thetai, lambda);   %plotting coordinates

A = zeros(n);   %container for absorption values

df = 25e-9;   %film thickness 25 nm

%transfer matrix method
k0 = @(omega) omega/c;   %free space wave number (1/m)

%full transfer matrix terms
M11 = @(t,omega,theta,epsf) cos(sqrt(epsf-(sin(theta))^2) * k0(omega) * t);
M12 = @(t,omega,theta,epsf) 1i*sin(sqrt(epsf-(sin(theta))^2)* k0(omega) * t)*sqrt((epsf-(sin(theta))^2)/(epsf)^2);
M21 = @(t,omega,theta,epsf) 1i*sin(sqrt(epsf-(sin(theta))^2) * k0(omega) * t)*sqrt(((epsf)^2)/(epsf-(sin(theta))^2));
M22 = @(t,omega,theta,epsf) cos(sqrt(epsf-(sin(theta))^2) * k0(omega) * t);
eta = @(theta) 1/cos(theta);

%reflection coefficient
r = @(t,omega,theta,epsf) (eta(theta)*M22(t,omega,theta,epsf) - eta(theta)*M11(t,omega,theta,epsf) - ((eta(theta)^2)*M12(t,omega,theta,epsf) - M21(t,omega,theta,epsf)))/(eta(theta)*M22(t,omega,theta,epsf) + eta(theta)*M11(t,omega,theta,epsf) - ((eta(theta)^2)*M12(t,omega,theta,epsf) + M21(t,omega,theta,epsf)));
%transmission coefficient
t = @(t,omega,theta,epsf) 2*eta(theta)/(eta(theta)*M22(t,omega,theta,epsf) + eta(theta)*M11(t,omega,theta,epsf) - ((eta(theta)^2)*M12(t,omega,theta,epsf) + M21(t,omega,theta,epsf)));

%compute absorption values
for i = 1:n
    for j = 1:n
        thetaij = deg2rad(T(i,j));   %current thetai
        lambdaij = L(i,j) * 1e-6;   %current lambda, convert μm to m
        omegaij = 2*pi*c/lambdaij;
        epsf = exp_eps(1,i);

        %reflection & transmission coefficients
        rij = r(df, omegaij, thetaij, epsf);
        tij = t(df, omegaij, thetaij, epsf);

        %absorption value assignment
        A(i,j) = 1 - abs(rij)^2 - abs(tij)^2;
    end
end

%plot
clf;
surf(T,L,A,'EdgeColor','None');   %3d plot of A
view(2);   %2d view
colormap hot
colorbar
xlim([t1,90])
ylim([0.250019531 1.684528687])
xlabel('Incident angle (degree)')
ylabel('Wavelength (μm)')
title('Absorption vs. incident angle and wavelength for df = 25 nm')

%max absorption
Amax = max(A,[],"all")
