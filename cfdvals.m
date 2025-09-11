%%==============%%

% Sameer Bajaj
% AS CFD Parameter Generator 2024-2025

% This program generates parameters needed for running rocket CFD from an
% imported openrocket flight CSV. After running the program, a table with
% all needed information will be copied to your clipboard (paste into 
% excel) and a graph of the initial layer height and cell count over time 
% will be generated

%%==============%%

clear
close all
clc

%%%% Set mode = "main" to run main program, "validate" to run equation validation

mode = "main";
switch mode
    case "main"
%% Import flight info
%{
Make sure your openrocket exports altitude (ft), vertical velocity (ft/s), 
mach number, reynolds number, temperature (F) & pressure (mbar). Also make
sure you turn flight events off. This is designed to work w/ default
imperical units in ORK.
%}
% Insert your file name below
flight_info = readmatrix("cfddata.csv");
% deletes rows after recovery
flight_info = flight_info(flight_info(:,2) >= 0, :);
% unit conversions
in2m = 0.0254;
ft2m = in2m*12;
mbar2Pa = 100;
flight_info(:, 1) = flight_info(:, 1) * ft2m;
flight_info(:, 2) = flight_info(:, 2) * ft2m;
flight_info(:, 5) = 5/9*(flight_info(:, 5) - 32) + 273.15;
flight_info(:, 6) = mbar2Pa * flight_info(:, 6);

altitude = flight_info(:, 1); %all rows, col 1
u_inf = flight_info(:,2); %flight velocity (m/s)
mach = flight_info(:,3);%mach number
Re_L = flight_info(:,4); %reynolds number
T_atm = flight_info(:,5); % (K)
P_atm = flight_info(:,6);% (Pa)

%% Rocket Properties
% Change these based on rocket geometry
fineness = 4; % Nosecone Fineness ratio
OD = 8.375; % (inches)
BT_length = 12.1; %boattail length
body_len = 135.42; % in
len = fineness*OD+BT_length+body_len;
len_m = len*in2m;

%% Fluid properties
% gamma is the c_p / c_v, it is called k in 105a
gamma = 1.4;
% Viscosity of air (sutherland's formula), modified below paper to work
% with Pa*s and degrees kelvin ⊆(⌒Ꮂ⌒)⊇
% https://www.grc.nasa.gov/www/k-12/airplane/viscosity.html
u_0 = 1.716*10^-5; % Pa*s
T_0 = 273.11; %K
S = 110.56;
mu = u_0 * ((T_atm/T_0).^1.5).*((T_0+S)./(T_atm+S));
rho = P_atm./(287.058 * T_atm);
%% Find C_f - skin friction coefficient
% (https://www.sjsu.edu/ae/docs/project-thesis/Ben.Hopkins-F22.pdf)
% skin roughness value (in)
k = 1.2*10^-3;
Re_cutoff = (mach >= 0.8) .* (44.62 * (len/k)^1.053 .* mach.^1.16) + ...
         (mach < 0.8)  .* (38.21 * (len/k)^1.053);
C_f = 0.455./(log10(min(Re_L, Re_cutoff)).^2.58.*(1+0.144.*mach.^2).^0.65);

% Find T_w - shear wall stress
tau_w = C_f.*(rho/2.*(u_inf.^2));

%% Find CFD input values

% y+: non-dimensional measure of what region of the boundary layer we are
% in. 
% 
% y+ << 1: you are right next to the wall and resolve the boundary layer 
% rigorously, don't use any wall treatment
%
% y+ < 5: you are in the viscous layer, and we approximate using a linear 
% equation. dont go past y+ = 5
%
% y+ < 30: in the buffer layer, less accurate but faster, use y+ = 25 if
% you get a low min orthogonal quality

y_plus = 5;
% r is your growth rate (usually 1.0-1.5)
r = 1.3;


% wall friction velocity
u_tau = sqrt(tau_w./rho);
% v is kinematic viscosity, equal to dynamic viscosity / density
v = mu./rho;
% y is first layer height (equation gives centroid, multiply by 2 to get
% height)
y = 2*y_plus*v./u_tau;
% N is the number of cells we need to use to approximate behavior, total
% boundary layer height = y*((1-r^N)/(1-r)).
height_BL = 0.382*len_m./(Re_L.^0.2);
N = log(1-(1-r).*height_BL./y)/log(r);

%% Plot Results
figure;
grid on; hold on;
yyaxis left
plot(altitude, y, 'Color', [39 116 174]/255, 'LineWidth', 1, ... 
    'LineStyle','-');
ylabel('Initial Layer Height (m)', 'Color', [39 116 174]/255)
ylim([0,prctile(y, 75)])
xlabel('Altitude (m)')
yyaxis right
plot(altitude, N, 'Color', [255 209 0]/255, 'LineWidth', 1, ... 
    'LineStyle', '-.');
ylabel('Number of Cells', 'Color', [0 0 0])
title('y and N vs Altitude')
legend('Layer Height', 'Cell Count')

%% Output values
[min_y_val, min_index] = min(y);
p25 = prctile(y, 25); 
p50 = prctile(y, 50);
[y_new, indices] = sort(y);
i25 = find(y_new >= p25, 1, 'first');
i50 = find(y_new >= p50, 1, 'first');
index_25 = indices(i25); index_50 = indices(i50);

output_table = {};
output_table{end + 1} = 'Rocket Geometry';
output_table{end + 1} = sprintf('OD\t%.1f\tin', OD);
output_table{end + 1} = sprintf('Length\t%.2f\tin', len);
output_table{end + 1} = sprintf('NC Fineness\t%.1f', fineness);
output_table{end + 1} = ' ';

output_table{end + 1} = 'Atmospheric Conditions (min, 25th, 50th y)';
output_table{end + 1} = sprintf('Altitude\t%.1f\t%.1f\t%.1f\tm', altitude(min_index), altitude(index_25), altitude(index_50));
output_table{end + 1} = sprintf('\t%.1f\t%.1f\t%.1f\tft', altitude(min_index)/ft2m, altitude(index_25)/ft2m, altitude(index_50)/ft2m);
output_table{end + 1} = sprintf('P_atm\t%.2E\t%.2E\t%.2E\tPa', P_atm(min_index), P_atm(index_25), P_atm(index_50));
output_table{end + 1} = sprintf('T_atm\t%.2f\t%.2f\t%.2f\tK', T_atm(min_index), T_atm(index_25), T_atm(index_50));
output_table{end + 1} = sprintf('mu\t%.2E\t%.2E\t%.2E\tPa*s', mu(min_index), mu(index_25), mu(index_50));
output_table{end + 1} = sprintf('rho\t%.2f\t%.2f\t%.2f\tkg/m^3', rho(min_index), rho(index_25), rho(index_50));
output_table{end + 1} = ' ';

output_table{end + 1} = 'Flow Properties (min, 25th, 50th y)';
output_table{end + 1} = sprintf('u_inf\t%.2f\t%.2f\t%.2f\tm/s', u_inf(min_index), u_inf(index_25), u_inf(index_50));
output_table{end + 1} = sprintf('Mach\t%.2f\t%.2f\t%.2f', mach(min_index), mach(index_25), mach(index_50));
output_table{end + 1} = sprintf('Re_L\t%.2E\t%.2E\t%.2E', Re_L(min_index), Re_L(index_25), Re_L(index_50));
output_table{end + 1} = sprintf('Re_cutoff\t%.2E\t%.2E\t%.2E', Re_cutoff(min_index), Re_cutoff(index_25), Re_cutoff(index_50));
output_table{end + 1} = sprintf('Cf\t%.8E\t%.8E\t%.8E', C_f(min_index), C_f(index_25), C_f(index_50));
output_table{end + 1} = sprintf('tau_w\t%.2f\t%.2f\t%.2f\tPa', tau_w(min_index), tau_w(index_25), tau_w(index_50));
output_table{end + 1} = sprintf('u_tau\t%.3f\t%.3f\t%.3f\tm/s', u_tau(min_index), u_tau(index_25), u_tau(index_50));
output_table{end + 1} = ' ';

output_table{end + 1} = 'CFD Input Values';
output_table{end + 1} = sprintf('Y+\t%.1f', y_plus);
output_table{end + 1} = sprintf('r\t%.2f', r);
output_table{end + 1} = sprintf('y1\t%.8E\t%.8E\t%.8E\tm', ...
    min_y_val, p25, p50);
output_table{end + 1} = sprintf('delta\t%.8E\t%.8E\t%.8E\tm', ...
    height_BL(min_index), height_BL(index_25), height_BL(index_50));
output_table{end + 1} = sprintf('N\t%d\t%d\t%d', ceil(N(min_index)), ceil(N(index_25)), ceil(N(index_50)));
output_table{end + 1} = ' ';

output_table{end + 1} = sprintf('NC Area\t%.7E\tm^2', pi*OD^2/4*in2m^2);
output_table{end + 1} = sprintf('Length\t%.2f\tm', len_m);
output_table{end + 1} = sprintf('C1: Ref Viscosity\t%.5E\tkg/m^3', u_0);
output_table{end + 1} = sprintf('C2: Ref Temp\t%.2f\tK', T_0);
output_table{end + 1} = sprintf('C3: Sutherland''s Constant\t%.2f', S);
output_table{end + 1} = sprintf('Velocity\t%.2f\t%.2f\t%.2f\tm/s', u_inf(min_index), u_inf(index_25), u_inf(index_50));
output_table{end + 1} = sprintf('Viscosity\t%.2E\t%.2E\t%.2E\tPa*s', mu(min_index), mu(index_25), mu(index_50));

final_output = strjoin(output_table, '\n');
clipboard('copy', final_output);
disp('Data has been copied to your clipboard.');





    case "validate"

%% Testing for coefficient of skin friction
Re_L = logspace(5,9, 40);
figure;
hold on; grid on;
set(gca, 'XScale', 'log')
ylim([0, 0.006])
xlim([10^5, 10^9]);
legend_out = {};
for mach = [0, .3, .7, .9, 1, 1.5, 2, 2.5, 3]
C_f = 0.455./(log10(Re_L).^2.58.*(1+0.144*mach^2).^0.65);
plot(Re_L, C_f, 'LineWidth', 1)
legend_out{end+1} = sprintf('Mach %.1f', mach);
end
xlabel('Re_L')
ylabel('Skin Friction Coefficient')
legend(legend_out);

%% Testing for cutoff reynolds number
lenk = logspace(2, 7, 50);
figure;
hold on; grid on;
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([10^5, 10^9])
ylim([10^2, 10^7])
legend_out = {};
for mach = [0, 1, 2, 3]
    Re_cutoff = (mach >= 0.8) .* (44.62 .* (lenk).^1.053 .* mach.^1.16) + ...
         (mach < 0.8)  .* (38.21 .* (lenk).^1.053);
plot(Re_cutoff, lenk, 'LineWidth', 1)
legend_out{end+1} = sprintf('Mach %.1f', mach);
end
xlabel('Re_{cutoff}')
ylabel('L / k')
legend(legend_out);
end