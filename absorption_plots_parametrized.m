%validation of absorption vs. conductivity and incident angle plot from liu et al.
%calculated using parametrized reflection/transmission coefficients

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
omega = 2*pi*f;   %angular frequency
w = 0.001;   %operating wavelength 1 mm

%constants
ep0 = 8.85419*10^-12;   %vacuum permittivity
mu0 = 4*pi*10^-7;   %vacuum permeability

%compute absorption values
for i = 1:n
    for j = 1:n
        t = deg2rad(T(i,j));   %current thetai
        s = S(i,j);   %current sigma

        beta = s*omega*mu0*df^2;
        gamma = sqrt(s/(omega*ep0))*cos(t);

        %reflection & transmission coefficients
        r = -(1i + gamma^2)*sin(sqrt(1i*beta)) / (2*sqrt(1i)*gamma*cos(sqrt(1i*beta)) + (-1i + gamma^2)*sin(sqrt(1i*beta)));
        t = -2*sqrt(1i)*gamma / (2*sqrt(1i)*gamma*cos(sqrt(1i*beta)) + (-1i + gamma^2)*sin(sqrt(1i*beta)));

        %absorption value assignment
        A(i,j) = 1 - abs(r)^2 - abs(t)^2;
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
