%DC conductivity of gold thin film

clear;

%retrieve experimental data
fileID = fopen('dcconductivity-vs-thickness.txt');
sdc_data = fscanf(fileID, '%f%f', [2 16]);
fclose(fileID);

%scale datapoints
for i = 1:16
    si = sdc_data(2,i);
    sdc_data(2,i) = si*10^8;
end

%fit to data (logarithmic)
sdc_fit = @(t) 0.0858*log(t) - 0.0708;   %conductivity (10^8 S/m) as a function of film thickness (nm)

n = 500;   %number of elements
df = linspace(20,180,n);   %film thickness (nm) -> converted to m later
sigma_dc = zeros([1 n]);

for i = 1:n
    dfi = df(1,i);
    sigma_dc(1,i) = sdc_fit(dfi)*10^8;
end

clf
hold on
scatter(sdc_data(1,:), sdc_data(2,:),'filled','b')   %experimental data
plot(df, sigma_dc,'b')   %fit to data
title('DC conductivity vs. film thickness')
xlabel('Film thickness (nm)')
ylabel('DC conductivity (S/m)')
hold off
legend('Measured','Fit to data',Location='northeastoutside')