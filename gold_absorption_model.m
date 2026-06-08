%analytical absorption model - gold

clear;

%number of elements
n = 500;

%independent variables
thetai = linspace(70, 90, n);   %incident angle (degrees), linear scale
df = linspace(10, 100, n);   %film thickness (nm), linear scale
lambda0 = linspace(600, 50000, n);   %operating free space wavelength (nm), linear scale

lambda1 = 900*10^-9;   %fixed wavelength for thickness sweep
omega1 = 2*pi*3*10^8/lambda1;   %corresponding omega
df1 = 1*10^-9;   %fixed thickness for wavelength sweep
%df1 = 1*10^-6;
theta1 = 89.8;   %fixed angle for thickness & wavelength sweep

[T1, Df] = meshgrid(thetai, df);   %plotting coordinates for varying thickness
[T2, L0] = meshgrid(thetai, lambda0);   %plotting coordinates for varying wavelength
[Df1, L01] = meshgrid(df, lambda0);   %plotting coordinates for varying thickness and wavelength
A1 = zeros(n); A2 = zeros(n); A3 = zeros(n);   %containers for absorption values

%Drude parameters
epsinf = 1.53;   %dielectric constant at infinite frequency
hbar = 6.582*10^-16;   %reduced planck's constant (eV*s)
omegap = 8.5/hbar;   %plasma frequency (1/s)

%damping rate
gammaep = 4.6*10^13;   %electron-phonon scattering (1/s)
gammas = 1.35*10^13;   %scattering on film surface (1/s)

R = 0.31;   %grain-boundary reflection coefficient
Lambda = 39*10^-9;   %mean free path for conducting electrons in bulk gold (m)
D = @(t) (1.522*10^-8)*log(t) + 2.847*10^-7;   %average crystalite size as a function of thickness (m)
alpha = @(t) Lambda*R/(D(t)*(1 - R));

%electron-grain boundary scattering rate (1/s)
gammagb = @(t) gammaep*(1/(1 - 1.5*alpha(t) + 3*alpha(t)^2 - 3*alpha(t)^3*log(1 + 1/alpha(t))) - 1);

gamma = @(t) gammaep + gammagb(t) + gammas;   %total damping parameter (1/s)

%dielectric function of film - Drude model with thickness correction
epsf = @(t,omega) epsinf - omegap^2/(omega^2 + gamma(t)^2) + 1i*gamma(t)*omegap^2/(omega*(omega^2 + gamma(t)));

%transfer matrix method
k0 = @(omega) omega/(3*10^8);   %free space wave number (1/m)

%simplified transfer matrix terms
% M11 = @(t,omega) cos(sqrt(epsf(t,omega)) * k0(omega) * t);
% M12 = @(t,omega) 1i*sin(sqrt(epsf(t,omega))* k0(omega) * t)/sqrt(epsf(t,omega));
% M21 = @(t,omega) 1i*sqrt(epsf(t,omega))*sin(sqrt(epsf(t,omega)) * k0(omega) * t);
% M22 = @(t,omega) cos(sqrt(epsf(t,omega)) * k0(omega) * t);

M11 = @(t,omega,theta) cos(sqrt(epsf(t,omega)-(sin(theta))^2) * k0(omega) * t);
M12 = @(t,omega,theta) 1i*sin(sqrt(epsf(t,omega)-(sin(theta))^2)* k0(omega) * t)*sqrt((epsf(t,omega)-(sin(theta))^2)/(epsf(t,omega))^2);
M21 = @(t,omega,theta) 1i*sin(sqrt(epsf(t,omega)-(sin(theta))^2) * k0(omega) * t)*sqrt(((epsf(t,omega))^2)/(epsf(t,omega)-(sin(theta))^2));
M22 = @(t,omega,theta) cos(sqrt(epsf(t,omega)-(sin(theta))^2) * k0(omega) * t);

eta = @(theta) 1/cos(theta);

%reflection coefficient
r = @(t,omega,theta) (eta(theta)*M22(t,omega,theta) - eta(theta)*M11(t,omega,theta) - ((eta(theta)^2)*M12(t,omega,theta) - M21(t,omega,theta)))/(eta(theta)*M22(t,omega,theta) + eta(theta)*M11(t,omega,theta) - ((eta(theta)^2)*M12(t,omega,theta) + M21(t,omega,theta)));
%transmission coefficient
t = @(t,omega,theta) 2*eta(theta)/(eta(theta)*M22(t,omega,theta) + eta(theta)*M11(t,omega,theta) - ((eta(theta)^2)*M12(t,omega,theta) + M21(t,omega,theta)));

% eps = zeros([1 n]);
% eps1 = zeros([1 n]);
% eps2 = zeros([1 n]);
% g = zeros([1 n]);
% Dt = zeros([1 n]);
% 
% for i = 1:n
%     %dfi = df(1,i)*10^-9;
%     dfi = 25*10^-9;
%     omegai = 2*pi*3*10^8 / (lambda0(1,i)*10^-9);
%     %omegai = 900*10^-9;
%     eps(1,i) = epsf(dfi,omegai);
%     eps1(1,i) = real(eps(1,i));
%     eps2(1,i) = imag(eps(1,i));
%     Di = D(dfi);
%     Dt(1,i) = Di;
% 
%     gi = gamma(dfi);
%     g(1,i) = gi;
% 
% end
% 
% %plot(lambda0,eps1)
% plot(lambda0,eps2)
% %xlim([200,2000])
% %ylim([-200,0])
% %plot(df, g)
% %plot(df,Dt)

%thickness sweep, fixed wavelength
for i = 1:n
    for j = 1:n
        thetaiij = deg2rad(T1(i,j));   %current thetai
        dfij = Df(i,j)*10^-9;   %current df

        rij = r(dfij,omega1,thetaiij);   %reflection coefficient
        tij = t(dfij,omega1,thetaiij);   %transmission coefficient
        Aij = 1 - abs(rij)^2 - abs(tij)^2;   %absorption
        A1(i,j) = Aij;   %assign absorption
    end
end

%plot absorption vs. film thickness and incident angle
figure;
surf(T2,Df,A1,'EdgeColor','None');   %3d plot of A1
view(2);   %2d view
colormap hot
colorbar
xlim([70,90])
xlabel('Incident angle (degree)')
ylabel('Film thickness (nm)')

%max absorption from thickness sweep
Amax_t = max(A1,[],"all")

%wavelength sweep, fixed thickness
for i = 1:n
    for j = 1:n
        thetaiij = deg2rad(T1(i,j));   %current thetai
        lambda0ij = L0(i,j)*10^-9;   %current lambda0
        omegaij = 2*pi*3*10^8/lambda0ij;   %corresponding omega

        rij = r(df1,omegaij,thetaiij);   %reflection coefficient
        tij = t(df1,omegaij,thetaiij);   %transmission coefficient
        Aij = 1 - abs(rij)^2 - abs(tij)^2;   %absorption
        A2(i,j) = Aij;   %assign absorption
    end
end

%plot absorption vs. wavelength and incident angle
figure;
surf(T1,L0,A2,'EdgeColor','None');   %3d plot of A2
view(2);   %2d view
colormap hot
colorbar
xlim([70,90])
xlabel('Incident angle (degree)')
ylabel('Operating wavelength (nm)')

%max absorption from wavelength sweep
Amax_w = max(A2,[],"all")


% %plot absorption vs. thickness and wavelength
% for i = 1:n
%     for j = 1:n
%         thetaiij = deg2rad(theta1);   %current thetai
%         dfij = Df1(i,j)*10^-9;   %current df
%         lambda0ij = L01(i,j)*10^-9;   %current lambda
%         omegaij = 2*pi*3*10^8/lambda0ij;   %corresponding omega
% 
%         rij = r(dfij,omega1,thetaiij);   %reflection coefficient
%         tij = t(dfij,omega1,thetaiij);   %transmission coefficient
%         Aij = 1 - abs(rij)^2 - abs(tij)^2;   %absorption
%         A3(i,j) = Aij;   %assign absorption
%     end
% end
% 
% %plot absorption vs. film thickness and incident angle
% figure;
% surf(Df1,L01,A3,'EdgeColor','None');   %3d plot of A1
% view(2);   %2d view
% colormap hot
% colorbar
% xlabel('Film thickness (nm)')
% ylabel('Operating wavelength (nm)')
% 
% %max absorption from thickness sweep
% Amax_tw = max(A3,[],"all")