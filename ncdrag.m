%{

Nosecone drag impulse evaluation
Sameer Bajaj

%}

% For ORK: Have your system output a csv with time, vertical velocity,
% mach number, temperature (F), pressure (mbar)

% For RASAERO: Export your aeroplot csv with the name "NC_len.csv", where
% NC_len is the length of the nosecone you are modelling

clear
close all
clc

% Input parameters for your nosecone:
num_layers = 4;
t_layer = 0.015;
OD = 6;
cm2in = 2.54;
in2m = 0.0254;
ft2m = in2m*12;
g = 9.81;
mbar2Pa = 100;
density_NC = 2.178*(100^3); %g/cm3 to g/m3
density_tip = 2.70 * 100^3; %g/cm3 to g/m3
N2lb = 1/4.4482216153;
tip_length = 1.5; %inches

mainloc = 'C:\Users\Owner\Downloads\NC_trials';
output_table = {};
header{1} = sprintf('Length (in)\tTotal Impulse (N*s)\tMass (lb)');
filenames = {dir(mainloc).name};
filenames = lower(filenames(:));
filenames = filenames(not(strcmp(filenames(:), '.') + ...
    strcmp(filenames(:), '..') + endsWith(filenames, '.cdx1')));
ORKfile = fullfile(mainloc, filenames{endsWith(filenames,'ork.csv')});
ORK_info = readmatrix(ORKfile);
ORK_info = ORK_info(ORK_info(:,3) > 0, :);
ORK_info(:, 2) = ORK_info(:, 2) * ft2m;
ORK_info(:, 3) = round(ORK_info(:, 3), 2);
ORK_info(:, 4) = 5/9*(ORK_info(:, 4) - 32) + 273.15;
ORK_info(:, 5) = mbar2Pa * ORK_info(:, 5);
time = ORK_info(:, 1);
v = ORK_info(:, 2);
mach = ORK_info(:, 3);
T_atm = ORK_info(:, 4);
P_atm = ORK_info(:, 5);
rho = P_atm./(287.058 * T_atm);
A = OD^2/4*pi*in2m^2; %m^2
mach_idx = round(mach*100);
RAS_outputs = filenames(~endsWith(filenames,'ork.csv'));
lengths = zeros(length(RAS_outputs), 1);
impulses = zeros(length(RAS_outputs), 1); 
for i = 1:length(RAS_outputs)
    NC_len = string(RAS_outputs{i});
    NC_len{1}(end-3:end) = [];
    NC_len = str2double(NC_len);
    % Calculate volume of nosecone
    vol_NC = getVol(NC_len, OD, t_layer*num_layers);
    t_tip = acos(1 - tip_length * 2 / NC_len);
    OD_tip = OD/sqrt(pi)*sqrt(t_tip-sin(2*t_tip)/2);
    vol_tip = getVol(tip_length, OD_tip, 0);
    mass_tip = vol_tip*in2m^3*density_tip/1000;
    disp(mass_tip)
    mass_nc = vol_NC*in2m^3*density_NC/1000; %convert from gram to kg
    mass = mass_tip + mass_nc;
    mass_force = mass*g;
    T = readtable(fullfile(mainloc, RAS_outputs{i}));
    strs = T.Var1;
    % get rid of all alphabetical elements
    mask = cellfun(@(s) all(~isstrprop(s, 'alpha'), 'all'), strs);
    strs = strs(mask);
    % turn string cell into a numerical cell
    C = cellfun(@(s) sscanf(s, '%f'), strs, 'UniformOutput', false);
    % get only arrays with more than 12 elements in C
    mask2 = cellfun(@(s) numel(s) >= 12, C);
    RAS_info = C(mask2);
    RAS_info = RAS_info(mach_idx);
    Cd = cellfun(@(s) s(4) - s(7), RAS_info);
    drag_force = 1/2 * A * Cd .* rho .* v .^ 2;
    F_tot = drag_force + mass_force;
    impulse = trapz(time, F_tot);
    output_table{end+1} = sprintf('%.2f\t%.5f\t%.4f', round(NC_len, 2), ...
        impulse, mass_force * N2lb);
    impulses(i) = impulse;
    lengths(i) = NC_len;
end
[lengths, idx_shift] = sort(lengths);
impulses = impulses(idx_shift);
output_table = output_table(idx_shift);
output_table = [header, output_table];
[min_impulse, min_impulse_idx] = min(impulses);
fprintf('Optimal NC Length Case: %.2f in\n', lengths(min_impulse_idx))
figure;
hold on
plot(lengths, impulses, 'Color', [39 116 174]/255, 'LineWidth', 1);
plot(lengths(min_impulse_idx), min_impulse, 'Marker', 'o', 'Color', [255 209 0]/255)
xlabel('Nosecone Length (in)');
ylabel('Total Impulse (N*s)');
title('Nosecone Drag Impulse Evaluation');
grid on;
final_output = strjoin(output_table, '\n');
clipboard('copy', final_output);
disp('Data has been copied to your clipboard.');
