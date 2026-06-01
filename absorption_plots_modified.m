%validation of absorption vs. conductivity and incident angle plot from liu et al.
%calculated using reflection/transmission coefficients in terms of transfer matrix elements

clear;

%number of elements
n = 1000;
%incident angle lower bound
t1 = 87;

%independent variables
thetai = linspace(t1, 90, n);   %incident angle (in degrees) on x-axis, linear scale
sigma = logspace(4, 7, n);   %conductivity on y-axis, logarithmic scale
%note - change endpoints of thetai and sigma intervals to get 'zoomed in' regions of plot (to directly match plots presented by liu et al.)
%change x-axis limits in line 52 to match thetai endpoints

[T,S] = meshgrid(thetai, sigma);   %plotting coordinates
A = zeros(n);   %container for absorption values

%film parameters
df = 10^-5;   %film thickness 10 μm
f = 0.3*10^12;   %operating freqeuncy 0.3 THz
Omega = 2*pi*f;   %angular frequency
w = 0.001;   %operating wavelength 1 mm

%constants
ep0 = 8.85419*10^-12;   %vacuum permittivity
mu0 = 4*pi*10^-7;   %vacuum permeability

%film dielectric
epsf = @(sig,omega) 1+1i*sig/(omega*ep0);

%transfer matrix method
k0 = @(omega) omega/(2*pi*3*10^8);   %free space wave number (1/m)

%simplified transfer matrix terms
M11 = @(t,sig,omega,theta) cos(sqrt(epsf(sig,omega)-(sin(theta))^2) * k0(omega) * t);
M12 = @(t,sig,omega,theta) 1i*sin(sqrt(epsf(sig,omega)-(sin(theta))^2)* k0(omega) * t)*sqrt((epsf(sig,omega)-(sin(theta))^2)/(epsf(sig,omega))^2);
M21 = @(t,sig,omega,theta) 1i*sin(sqrt(epsf(sig,omega)-(sin(theta))^2) * k0(omega) * t)*sqrt(((epsf(sig,omega))^2)/(epsf(sig,omega)-(sin(theta))^2));
M22 = @(t,sig,omega,theta) cos(sqrt(epsf(sig,omega)-(sin(theta))^2) * k0(omega) * t);

eta = @(theta) 1/cos(theta);

%reflection coefficient
r = @(t,sig,omega,theta) (eta(theta)*M22(t,sig,omega,theta) - eta(theta)*M11(t,sig,omega,theta) - ((eta(theta)^2)*M12(t,sig,omega,theta) - M21(t,sig,omega,theta)))/(eta(theta)*M22(t,sig,omega,theta) + eta(theta)*M11(t,sig,omega,theta) - ((eta(theta)^2)*M12(t,sig,omega,theta) + M21(t,sig,omega,theta)));
%transmission coefficient
t = @(t,sig,omega,theta) 2*eta(theta)/(eta(theta)*M22(t,sig,omega,theta) + eta(theta)*M11(t,sig,omega,theta) - ((eta(theta)^2)*M12(t,sig,omega,theta) + M21(t,sig,omega,theta)));


%compute absorption values
for i = 1:n
    for j = 1:n
        th = deg2rad(T(i,j));   %current thetai
        s = S(i,j);   %current sigma

        %reflection & transmission coefficients
        rij = r(df,s,Omega,th);
        tij = t(df,s,Omega,th);

        %absorption value assignment
        A(i,j) = 1 - abs(rij)^2 - abs(tij)^2;
    end
end

%plot
figure;
surf(T,S,A,'EdgeColor','None');   %3d plot of A
view(2);   %2d view
yscale log
colormap hot
colorbar
xlim([t1,90])
xlabel('Incident angle (degree)')
ylabel('Conductivity (S/m)')

%max absorption
Amax = max(A,[],"all")
