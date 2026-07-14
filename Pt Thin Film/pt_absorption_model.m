%analytical absorption model - platinum

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
t1 = 60;   %incident angle lower bound
thetai = linspace(t1, 90, n);   %incident angle (in degrees) on x-axis, linear scale
lambda = linspace(0.19, 3.5, n);   %wavelengths 0.19 to 3.5 μm -> converted to m later
[T,L] = meshgrid(thetai, lambda);   %plotting coordinates

%containers
Epsf = zeros([1 n]);   %dielectric function of film

for i = 1:n
    lambdai = lambda(1,i) * 1e-6;   %convert μm to m
    omega = 2*pi*c/lambdai;
    
    %Lorentz oscillators
    lo = 0;
    for j = 1:4
        lo = lo + Omegaj(1,j)^2 / (omegaj(1,j)^2 - omega^2 - 1i*gammaj(1,j)*omega);
    end

    %Drude-Lorentz
    Epsf(1,i) = epsinf - Omegap^2 / (omega^2 + 1i*Gamma*omega) + lo;
end


%absorption calculation
A = zeros(n);   %container for absorption values
df = 80e-9;   %film thickness 80 nm

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
        epsf = Epsf(1,i);

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
ylim([0.19 3.5])
xlabel('Incident angle (degree)')
ylabel('Wavelength (μm)')
title('Absorption vs. incident angle and wavelength for df = 80 nm')

%max absorption
Amax = max(A,[],"all")
