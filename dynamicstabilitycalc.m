%{

Dynamic Stability Calculations
By: Sameer Bajaj

%}


%% Instructions for use!!!
% GENERAL:
% create a folder, copy its path and put it in below as 'mainloc'
% create folders within that folder for each case you want to test, can
% name them whatever you want.
% export all csvs as directed below into your folders.
% make sure while running the program you don't have the folders open.

% IN OPENROCKET:
% This script is designed to work with default imperical units
% Export time (s), altitude (ft), vertical velocity (ft/s), motor mass (oz)
% longitudal moment of inertia (lb-ft^2), 
% rotational moment of inertia (lb-ft^2), CG, location (in), stability
% margin, mach, reference area (in^2), air temp (F), air press (mbar)

% IF!! you want to get results including body tubes, go to tools -->
% component analysis and write down the CP locations of body tubes in
% inches in a txt or csv file called "BT". This doesn't change w mach so dw
% abt that.

% IN RASAERO:
% Take your fullrocket file, and create subpart files by removing
% components below the nosecone. Export your aeroplot results as a csv. 
% the aeroplot csv with just the nosecone should be called "RAS_1", the 
% one w nosecone and fins is "RAS_2" and the one with everything can be 
% "RAS_full" or "RAS_x", where x is the last number. This script is also
% designed to work if you wanna include body tubes. If you do, just have
% the file with NC and highest body tube be "RAS_2", and so on.

clear
close all
clc
output_table = {};
% Conversions
in2m = 0.0254;
ft2m = in2m*12;
mbar2Pa = 100;
lb2kg = 0.45359237;
oz2lb = 1/16;

% Info
L_ne = 173*in2m; % length of nozzle exit from NC
rail_len = 40*ft2m; % length of the launch rail
% enter your main folder location
mainloc = 'C:\Users\Owner\Documents\MATLAB\stab_calc_files';
% get the names of folders in the directory.
foldernames = {dir(mainloc).name};
% removes the beginning two folder names which are just . and ..
foldernames = foldernames(not(strcmp(foldernames(:), '.') + ...
    strcmp(foldernames(:), '..')));

output_table{end + 1} = sprintf('\tOTR damping\tOTR stability\t\tMin\t25th\tMedian\t75th');
for i = 1:length(foldernames)
    folder = foldernames{i};
    filenames = {dir(fullfile(mainloc, folder)).name};
    % removes the beginning two file names which are just . and ..
    filenames = filenames(not(strcmp(filenames(:), '.') + ...
        strcmp(filenames(:), '..')));
    filenames = lower(filenames(:));
    ORKfile = fullfile(mainloc, folder, filenames{endsWith(filenames,'ork.csv')});
    % fullrocketfile = fullfile(folder, filenames{endsWith(filenames, 'full.csv')});
    ORK_info = readmatrix(ORKfile);
    % full_RAS_info = readmatrix(fullrocketfile);
    ORK_info = ORK_info(ORK_info(:,3) > 0, :);
    ORK_info = ORK_info(~isnan(ORK_info(:,8)), :);
    ORK_info(:, 2) = ORK_info(:, 2)*ft2m;
    ORK_info(:, 3) = ORK_info(:, 3)*ft2m;
    ORK_info(:, 4) = ORK_info(:, 4)*oz2lb*lb2kg;
    ORK_info(:, 5) = ORK_info(:, 5)*ft2m^2*lb2kg;
    ORK_info(:, 6) = ORK_info(:, 6)*ft2m^2*lb2kg;
    ORK_info(:, 7) = ORK_info(:, 7) * in2m;
    ORK_info(:, 9) = round(ORK_info(:, 9), 2);
    ORK_info(:, 10) = ORK_info(:, 10)*in2m^2;
    % convert temperature from farenheit to kelvin
    ORK_info(:, 11) = 5/9*(ORK_info(:, 11) - 32) + 273.15;
    ORK_info(:, 12) = mbar2Pa * ORK_info(:, 12);
    % load in variables from 
    time = ORK_info(:,1); %s
    altitude = ORK_info(:,2); %m
    v = ORK_info(:, 3); %m/s
    m_dot = abs(diff(ORK_info(:, 4))./diff(time)); %kg/s
    m_dot(end+1) = 0;
    I_L = ORK_info(:,5); %kg*m^2
    I_r = ORK_info(:,6); %kg*m^2
    CG = ORK_info(:,7); % in
    static_stab = ORK_info(:,8);
    mach = ORK_info(:, 9);
    A_ref = ORK_info(:, 10);
    T_atm = ORK_info(:,11);
    P_atm = ORK_info(:,12);
    % gets our ras aero indices from the rounded mach numbers in openrocket
    mach_indices = round(mach*100);     
    rho = P_atm./(287.058 * T_atm);
    
    % Calculate corrective moment coefficient
    %% Solve rocket parts
    % get all files that start with "ras"
    AS_files = filenames(startsWith(filenames, 'ras'));
    AS_file_init = fullfile(mainloc, folder, AS_files{1});
    init_AS_data = readmatrix(AS_file_init);
    init_AS_data = init_AS_data(mach_indices, :);
    % initialize CP val and CNa val vectors to have the same rows as our
    % flight profile, and a column for every ras aero file csv
    CP_vals = zeros(length(mach_indices), length(AS_files));
    CNa_vals = CP_vals;
    % convert CP vals from inches to meters
    init_AS_data(:, 13) = init_AS_data(:, 13) * in2m;
    % first row of CP vals will be our first set of data
    CP_vals(:, 1) = init_AS_data(:, 13);
    CNa_vals(:, 1) = init_AS_data(:, 12);
    for j = 2:length(AS_files)
        cur_AS_file = fullfile(mainloc, folder, AS_files{j});
        cur_AS_data = readmatrix(cur_AS_file);
        cur_AS_data = cur_AS_data(mach_indices, :);
        cur_AS_data(:, 13) = cur_AS_data(:, 13) * in2m;
        CNa_vals(:, j) = cur_AS_data(:, 12) - sum(CNa_vals, 2);
        if length(AS_files) > 3 && j > 1 && j < length(AS_files) - 1
            BT_file = fullfile(mainloc, folder, filenames{startsWith(filenames, 'bt')});
            BT_data = readmatrix(BT_file) * in2m;
            cur_CP = BT_data(j - 1);
        else
            cur_CP = (cur_AS_data(:, 13).*cur_AS_data(:, 12) - ... 
                sum(CNa_vals.*CP_vals, 2))./CNa_vals(:, j);
            % ^ based off barrowman equations
            % z_tot*cna*tot = z_i*cna_i + sum(z_k*cna_k) from k = 1 to k = i -1
            cur_CP(isnan(cur_CP)) = cur_CP(find(~isnan(cur_CP), 1));
        end
        CP_vals(:, j) = cur_CP;
    end
    CP_full = cur_AS_data(:, 13); %m
    CNa_full = cur_AS_data(:, 12); %(1/rad)
    C_1 = rho/2.*v.^2.*A_ref.*CNa_full.*(CP_full-CG);
    C_2a = rho/2.*v.*A_ref.*sum(CNa_vals.*(CP_vals-CG).^2, 2);
    C_2r = m_dot.*(L_ne-CG).^2;
    C_2 = C_2r + C_2a;
    damp_rat = C_2./(2*sqrt(C_1.*(I_L+I_r)));
    OTR_damp_rat = damp_rat(find(altitude > rail_len, 1));
    OTR_static_stab = static_stab(find(altitude > rail_len, 1));
    [burnout_idx, trash] = find(m_dot == 0, 1);
    fig = figure;
    grid on; hold on;
    plot(altitude/ft2m, damp_rat, 'Color', [39 116 174]/255, 'LineWidth', 1);
    % starts from the rail
    xlim([rail_len/ft2m inf])
    xline(altitude(burnout_idx)/ft2m, '-', 'Burnout',  ... 
        'LabelVerticalAlignment', 'bottom', 'LineWidth', 1, 'Color', [200 0 0]/255)
    xlabel('Altitude (ft)')
    ylabel('Damping Ratio')
    title(sprintf('%s Damping Ratio', string(folder)))
    img_file = fullfile(mainloc, folder, sprintf('%s_damping_ratio.jpg', string(folder)));
    exportgraphics(fig, img_file);
    output_table{end + 1} = sprintf('%s\t%.5f\t%.5f\t\t%.5f\t%.5f\t%.5f\t%.5f', string(folder), OTR_damp_rat, OTR_static_stab, min(damp_rat), prctile(damp_rat, 25), median(damp_rat), prctile(damp_rat, 75));
end
final_output = strjoin(output_table, '\n');
clipboard('copy', final_output);
disp('Data has been copied to your clipboard.');