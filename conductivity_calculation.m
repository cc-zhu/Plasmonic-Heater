%computing conductivity from dielectric function (model verification)
%calculating values for each element in the loop instead of defining functions

clear;

%number of elements
n = 500;

%independent variables
% thetai = linspace(0, 90, n);   %incident angle (degrees), linear scale
df = linspace(20, 180, n);   %film thickness (nm), linear scale
lambda0 = linspace(200, 2000, n);   %operating free space wavelength (nm), linear scale

lambda1 = 1*10^-6;   %fixed wavelength for thickness sweep
omega1 = 2*pi*3*10^8/lambda1;   %corresponding omega
df1 = 40*10^-9;   %fixed thickness for wavelength sweep
%df1 = 1*10^-6;

% [T1, Df] = meshgrid(thetai, df);   %plotting coordinates for varying thickness and incident angle
% [T2, L0] = meshgrid(thetai, lambda0);   %plotting coordinates for varying wavelength and incident angle
% A1 = zeros(n); A2 = zeros(n);   %containers for absorption values


%constant definitions
ep0 = 8.854 * 10^-12;   %vacuum permittivity
c = 3*10^8;   %speed of light

%Drude parameters
epsinf = 1.53;   %dielectric constant at infinite frequency
hbar = 6.582*10^-16;   %reduced planck's constant (eV*s)
omegap = 8.5/hbar;   %plasma frequency (1/s)

%damping rate
gammaep = 4.6*10^13;   %electron-phonon scattering (1/s)
gammas = 1.35*10^13;   %scattering on film surface (1/s)

R = 0.31;   %grain-boundary reflection coefficient
Lambda = 39*10^-9;   %mean free path for conducting electrons in bulk gold (m)

%containers
eps = zeros([1 n]);   %dielectric function values
eps1 = zeros([1 n]);   %real part of dielectric function
eps2 = zeros([1 n]);   %imaginary part of dielectric function
g = zeros([1 n]);   %damping parameter values
Dt = zeros([1 n]);   %average crystallite size values

%directly rearrange dielectric for conductivity
sigmar1 = zeros([1 n]);   %real part of conductivity
simgai1 = zeros([1 n]);   %imaginary part of conductivity

%complex conductivity calculation with separated components
sigmar2 = zeros([1 n]);
simgai2 = zeros([1 n]);

%dielectric & conductivity vs thickness
for i = 1:n
    t = df(i)*10^-9;   %current thickness (nm converted to m)
    omega = omega1;   %fixed frequency
    
    D = (1.522*10^-8)*log(t) + 2.847*10^-7;   %average crystallite size (m)
    Dt(i) = D;

    alpha = Lambda*R/(D*(1 - R));
    gammagb = gammaep*(1/(1 - 1.5*alpha + 3*alpha^2 - 3*(alpha^3)*log(1 + 1/alpha)) - 1);   %electron-grain boundary scattering rate (1/s)
    gamma = gammaep + gammagb + gammas;   %total damping parameter (1/s)
    g(i) = gamma;

    epsf = epsinf - omegap^2/(omega^2 + 1i*gamma*omega);   %dielectric function
    eps(i) = epsf;   %assign dielectric
    eps1(i) = real(eps(i));   %assign real part
    eps2(i) = imag(eps(i));   %assign imaginary part
    
    %directly rearrange for sigma
    sigma = -1i*omega*ep0*(epsf - 1);
    sigmar1(i) = real(sigma);   %real conductivity
    sigmai1(i) = imag(sigma);   %imaginary conductivity

    %solve for real and imaginary parts separately
    % sigmar2(i) = omega*ep0*eps2(i);
    % sigmai2(i) = -omega*ep0*(eps1(i)-1);
end

%plot(df,eps2)
plot(df,sigmar1)
%plot(df,sigmai2)

%dielectric & conductivity vs wavelength
for i = 1:n
    lambda = lambda0(i)*10^-9;   %current wavelength (nm converted to m)
    omega = 2*pi*c/lambda;   %current frequency
    t = df1;   %fixed thickness
    
    D = (1.522*10^-8)*log(t) + 2.847*10^-7;   %average crystallite size (m)
    Dt(i) = D;

    alpha = Lambda*R/(D*(1 - R));
    gammagb = gammaep*(1/(1 - 1.5*alpha + 3*alpha^2 - 3*(alpha^3)*log(1 + 1/alpha)) - 1);   %electron-grain boundary scattering rate (1/s)
    gamma = gammaep + gammagb + gammas;   %total damping parameter (1/s)
    g(i) = gamma;

    epsf = epsinf - omegap^2/(omega^2 + 1i*gamma*omega);   %dielectric function
    eps(i) = epsf;   %assign dielectric
    eps1(i) = real(eps(i));   %assign real part
    eps2(i) = imag(eps(i));   %assign imaginary part

    

    % sigmar(i) = real(sigma);   %real conductivity
    % sigmai(i) = imag(sigma);   %imaginary conductivity
end

%plot(lambda0,eps2)


%previous functions
% D = @(t) (1.522*10^-8)*log(t) + 2.847*10^-7;
% alpha = @(t) Lambda*R/(D*(1 - R));
% gammagb = @(t) gammaep*(1/(1 - 1.5*alpha(t) + 3*alpha(t)^2 - 3*alpha(t)^3*log(1 + 1/alpha(t))) - 1);
% gamma = @(t) gammaep + gammagb(t) + gammas;
% epsf = @(t,omega) epsinf - omegap^2/(omega^2 + gamma(t)^2) + 1i*gamma(t)*omegap^2/(omega*(omega^2 + gamma(t)));
