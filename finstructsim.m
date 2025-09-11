   %% Clear Cache
clear all %#ok<*call>
close all 
clc
%% Set parameters for output
to_plot = true;
FinSim = false;
%% Set fin coords & thicknesses

rc = 11;
tc = 4;
semispan = 6.5;
sweep = 7;

M = [
0	0;
sweep	semispan;
sweep+tc	semispan;
rc	0;
];
num_layers = 4;
t_skin = 0.028*2*num_layers; t_core = .35; thick = (t_skin+t_core); 
MACHNUMBER = 0.5; body_dia = 0;
x_1 = M(1,1); x_2 = M(2,1); x_3 = M(3,1); x_4 = M(4,1);
y_1 = M(1,2); y_2 = M(2,2); y_3 = M(3,2); y_4 = M(4,2);

%% Shear Moduli
% define material properties
G_f = 231e9/2/1.5; G_m = 2e9;
E_f = 231e9; E_m = 6e9;
G_1 = 2.47e9/(2*1.38); G_2 = 1.04973e+8/(2*1.25);
E_1 = 2.47e9; E_2 = 1.04973e+8; 

% nylon and foam density (lb/in3)
p_1 = 0.036488571129231; p_2 = 0.00271; p_skin = 0.057803677036405;

% atmospheric pressure
P0 = 101325;

%% Find Fin Parameters
semispan = max(M(:, 2)); tc = sqrt((x_3 - x_2)^2 + (y_3 - y_2)^2); 
rc = x_4; 
x = M(:, 1); y = M(:, 2);
% centroid calculations
fin = polyshape(x, y);
[Cx, Cy] = centroid(fin);
area = area(fin);
% Quarterline
x0 = 0.25*x_4; y0 = 0.25*y_4; x1 = x_2 + 0.25*(x_3 - x_2);
y1 = y_2 + 0.25*(y_3 - y_2); slope = ((y1-y0)/(x1-x0));
distance = Cx-((Cy-y0)/slope + x0);
epscillon = distance/rc;
Y = 24*epscillon*1.4/pi*P0;
AR = semispan^2/area;
lamb = tc/rc;

%% Plot airfoil shape
rc2 = rc-t_skin; tc2 = tc - t_skin;
[x_skin_rc, AF_skin_rc] = AF_plot(thick, rc, 4000);
[x_core_rc, AF_core_rc] = AF_plot(t_core, rc2, 4000);
x_core_rc = x_core_rc + t_skin/2;
[x_core_tc, AF_core_tc] = AF_plot(t_core,tc2,4000);
[x_skin_tc, AF_skin_tc] = AF_plot(thick,tc,4000);

%% Find mass values
avg_core_area = (trapz(x_core_tc, AF_core_tc) + trapz(x_core_rc, AF_core_rc));
avg_skin_area = (trapz(x_skin_rc, AF_skin_rc)-trapz(x_core_rc, AF_core_rc)) + ...
    (trapz(x_skin_tc, AF_skin_tc)-trapz(x_core_tc, AF_core_tc));
V_core = avg_core_area * semispan;
V_skin = avg_skin_area * semispan;
mass_mat1 = V_core * p_1 + V_skin * p_skin;
mass_mat2 = V_core * p_2 + V_skin * p_skin;
% average thickness across airfoils
area_core = 2*trapz(x_core_rc, AF_core_rc);
avg_core = area_core/rc2;

%% Find equivalent shear / youngs moduli
E_skin = HalpinTsai(E_m, E_f); G_skin = HalpinTsai(G_m, G_f);
G_fin1 = (thick)/(t_skin/G_skin+avg_core/G_1); G_fin2 = (thick)/(t_skin/G_skin+avg_core/G_2);
% G_fin1 = G_1*V_core/(V_skin+V_core)+G_skin*t_skin/thick; G_fin2 = G_2*V_core/(V_skin+V_core)+G_skin*t_skin/thick;
E_fin1 = Equilize(E_m, E_f, E_1, V_core, V_skin, 0); E_fin2 = Equilize(E_m, E_f, E_2, V_core, V_skin, 0);
G_fin1 = Equilize(G_m, G_f, G_1, V_core, V_skin, 0); G_fin2 = Equilize(G_m, G_f, G_2, V_core, V_skin, 0);
% E_fin1 = (E_skin*t_skin+E_1*t_core)/thick; E_fin2 = (E_skin*t_skin + E_2 * t_core)/thick;
%% Plot mach vs height
airp = readmatrix('rocketair.csv');
flut1 = sqrt(G_fin1/((Y*AR^3)/((thick/rc)^3*(AR+1.65))*((lamb+1)/2)));
flut2 = sqrt(G_fin2/((Y*AR^3)/((thick/rc)^3*(AR+1.65))*((lamb+1)/2)));
Pressure = 10^2*airp(:,4);
Temp = airp(:,3);
P1 = Pressure(1)./Pressure; 
sound = airp(:,5);
velocity = airp(:,2);
mach1 = flut1.*sqrt(P1);
mach2 = flut2.*sqrt(P1);
v_f1 = sound.*mach1;
v_f2 = sound.*mach2;
h = airp(:,1);


%% Find altitude indices + volumes
max_height = 15000; altitudes = [];
for i = 0:3000:max_height
    altitudes = [altitudes, find(h == i)];
end

%% Divergence Calculations!

sweang = atan((sweep+tc/2-rc/2)/semispan);
% Get airfoil equations (corrected for effective root)
% Find effective root chord

slope_mid = tan(-sweang);
if x_3 == rc
    fin_type = 2;
end
if sweep ~= 0
    slope_mid = tan(-sweang);
    x_1_rc = slope_mid*rc/2/(slope_mid-y_2/x_2);
    y_1_rc = y_2/x_2*x_1_rc;
    x_2_rc = (slope_mid*rc/2-y_3/(x_3-x_4)*rc)/(slope_mid-y_3/(x_3-x_4));
    y_2_rc = slope_mid*(x_2_rc-rc/2);
    if fin_type == 2
        x_2_rc = rc;
        y_2_rc = slope_mid*(x_2_rc-rc/2);
    end
    rc_eff = sqrt((y_2_rc-y_1_rc)^2+(x_2_rc-x_1_rc)^2);
    
    
    % Find effective tip chord
    
    x_1_tc = (slope_mid*(sweep+tc/2)-semispan)/(slope_mid-y_2/x_2);
    y_1_tc = y_2/x_2*x_1_tc;
    x_2_tc = (-y_3/(x_3-x_4)*rc-semispan+slope_mid*(sweep+tc/2))/(slope_mid-y_3/(x_3-x_4));
    y_2_tc = slope_mid*(x_2_tc-sweep-tc/2)+semispan;
    if fin_type == 2
        x_2_tc = tc+sweep;
        y_2_tc = slope_mid*(x_2_tc-sweep-tc/2)+semispan;
    end
    tc_eff = sqrt((y_2_tc-y_1_tc)^2+(x_2_tc-x_1_tc)^2);
    lamb = tc_eff/rc_eff;
    
    % Get equation for airfoil
    [x_eff, AF_eff] = AF_plot(thick, rc_eff, 4000);
else
    x_eff = x_skin_rc; AF_eff = AF_skin_rc; rc_eff = rc; y_2_rc = -abs(rc*sin(sweang))/2;
end
% Prandtl–Glauert correction
m_o = (2*pi/sqrt(1-MACHNUMBER^2));
% Aspect ratio for this is terms of span^2 / total area
AR = (2*semispan)^2/(area);
edge = sqrt(1+4*(cos(sweang))^2/(AR^2));
m_e = m_o*cos(sweang)/edge;
a_airfoil = 2*trapz(x_eff, AF_eff)/(12^2);
% calculating ratio of segment to thickness
s = 2*arclength(x_eff, AF_eff, 'spline')/12;
% Uses MIT torsion approximation for thin open bodies
c_avg = area/semispan; tau = thick/rc_eff;
I = 8/3*trapz(x_eff, AF_eff.^3);
K_I = 1/(rc_eff^4*tau^3)*I;
I1 = K_I * c_avg^4*16/((1+lamb)^4)*tau^3;
G_fin1 = G_fin1 / 6895; G_fin2 = G_fin2 / 6895;
e_1 = .25; L = sqrt(semispan^2+(x_2+tc*0.5-0.5*rc)^2);

% lb per sq in
E_fin1 = E_fin1 / 6895; E_fin2 = E_fin2  / 6895;
% Sets torsional constant equal to moment of inertia for the airfoil
J = I;
GJ1 = G_fin1 * J; EI1 = E_fin1 * I;
GJ2 = G_fin2 * J; EI2 = E_fin2 * I;
% Find constant Ks

if lamb >= 0.2 && lamb <=0.5
    K_1 = 2.81-0.13*(lamb-0.2)/0.3;
    K_2 = 0.614-0.117*(lamb-0.2)/0.3;
elseif lamb > 0.5 && lamb <= 1
    K_1 = 2.74-0.27*(lamb-0.5)/0.5;
    K_2 = .497-0.107*(lamb-0.5)/0.5;
elseif lamb > 1 && lamb <= 1.5
    K_1 = 2.47 - 0.23*(lamb-1)/0.5;
    K_2  = 0.39 - 0.064*(lamb-1)/0.5;
else
    error('Taper ratio out of acceptable range!')
end


q_d1 = (GJ1)/(m_e*rc_eff*L^3*(cos(sweang))^2)*(L/e_1/rc_eff)*(K_1/(1+K_2*(GJ1)/(EI1)*L/(e_1*rc_eff)*tan(sweang)));
q_d2 = (GJ2)/(m_e*rc_eff*L^3*(cos(sweang))^2)*(L/e_1/rc_eff)*(K_1/(1+K_2*(GJ2)/(EI2)*L/(e_1*rc_eff)*tan(sweang)));
q_d1 = q_d1 * (12^2); q_d2 = q_d2 * (12^2);
rho = 1/515.4 * Pressure./(287 .* (Temp + 273.15));
v_d1 = sqrt(q_d1./(rho./2));
v_d2 = sqrt(q_d2./(rho./2));

rc = x_4;
if to_plot
    figure(1);
    subplot(3,1,[2 3]);
    plot (h, mach1, 'r', 'DisplayName', 'MAT 1 Flutter Speed')
    hold on;
    plot(h, mach2, 'b', 'DisplayName', 'MAT 2 Flutter Speed')
    plot(h, velocity./sound, 'g', 'DisplayName', 'Flight Speed', 'LineWidth', 1)
    xlabel('Altitude (ft)')
    ylabel('Flutter speed (mach)')
    title('Fin Flutter Speed')
    xlim([0 inf])
    legend;
    grid on;
    % Plot the fin shape
    subplot(3,2,1)
    fill(M(:,1), M(:,2), 'w')
    title('Fin Shape')
    hold on;
    plot(Cx, Cy, 'ro');
    line([x0 x1], [y0, y1])
    hold on;
    if sweep~=0
    line([x_1_rc x_2_rc], [y_1_rc y_2_rc])
    line([x_1_tc x_2_tc], [y_1_tc y_2_tc])
    end
    line([rc/2, sweep+tc/2], [0, semispan])
    subplot(3,2,2)
    plot(x_core_rc, AF_core_rc, 'b')
    hold on;
    plot(x_core_rc, -AF_core_rc, 'b') 
    plot(x_skin_rc, AF_skin_rc, 'r')
    plot(x_skin_rc, -AF_skin_rc, 'r')
    title('Airfoil')
    axis equal;
    figure(3);
    plot (h, v_d1, 'r', 'DisplayName', 'MAT 1 Divergence Speed')
    hold on;
    plot(h, v_d2, 'b', 'DisplayName', 'MAT 2 Divergence Speed')
    plot(h, velocity, 'g', 'DisplayName', 'Flight Speed', 'LineWidth', 1)
    xlabel('Altitude (ft)')
    ylabel('Divergence speed (ft/s)')
    title('Fin Divergence Speed')
    xlim([0 inf])
    legend;
    grid on;
end

%% Factor of safety calculations
FOS = 1.6;
MOS1_d = min(abs(v_d1./FOS./velocity))-1;
MOS1_f = min(abs(v_f1./FOS./velocity))-1;
MOS2_f = min(abs(v_f2./FOS./velocity))-1;
MOS2_d = min(abs(v_d2./FOS./velocity))-1;
o_1 = E_fin1/(2*G_fin1) - 1; o_2 = E_fin2/(2*G_fin2) - 1;
%% print pastable results
if FinSim
    mat1_data = {}; mat2_data = {};
    mat1_data{end+1} = sprintf('%s\t%s\t%s', 'Altitude', 'Flutter', 'Divergence');
    for k = 1:length(altitudes)
        mat1_data{end + 1} = sprintf('%d\t%.1f\t%.1f', abs(h(altitudes(k))), v_f1(altitudes(k)), v_d1(altitudes(k)));
    end

    final_output = strjoin(mat1_data, '\n');
    clipboard('copy', final_output);
    disp('Data for material 1 has been copied to your clipboard!');

else
    compare_data = {};
    compare_data{end+1} = sprintf('\t%s\t\t%s', 'MAT1', 'MAT2');
    compare_data{end+1} = sprintf('%s\t%s\t%s\t%s\t%s\t', 'Altitude', 'Flutter', 'Divergence', ... 
        'Flutter', 'Divergence');
    for k = 1:length(altitudes)
        compare_data{end + 1} = sprintf('%d\t%.1f\t%.1f\t%.1f\t%.1f', abs(h(altitudes(k))), v_f1(altitudes(k)), v_d1(altitudes(k)), ...
            v_f2(altitudes(k)), v_d2(altitudes(k)));
    end
    compare_data{end + 1} = sprintf('\n%s\t%.2f\t\t%.2f', 'Mass (lb)', mass_mat1, mass_mat2);
    compare_data{end + 1} = sprintf('%s\t%.2f\t%.2f\t%.2f\t%.2f', 'MOS', MOS1_f, MOS1_d, ...
        MOS2_f, MOS2_d);
    final_output = strjoin(compare_data, '\n');
    clipboard('copy', final_output);
    disp('Data for both materials have been copied to your clipboard!');
end
























%% Functions

function N_skin = HalpinTsai(N_m, N_f)
    eps = 1.5; Vf = 0.50;
    n = (N_f/N_m - 1)/((N_f/N_m)+ eps);
    N_skin = N_m * (1+eps*n*Vf)/(1-n*Vf);
end
function N_fin = Equilize(N_m, N_f, N_core, V_core, V_skin, N_skin)
    if N_skin == 0
        eps = 1.5; Vf = 0.50;
        n = (N_f/N_m - 1)/((N_f/N_m)+ eps);
        N_skin = N_m * (1+eps*n*Vf)/(1-n*Vf);
    end    
    N_fin = (N_core*V_core+N_skin*V_skin)/(V_skin+V_core);
    
end
function [x, AF_plot] = AF_plot(thickness, chord, accuracy)
    x = linspace(0, chord, accuracy);
    AF_plot = 5 * (thickness) * (0.2969 * sqrt(x / chord) - 0.126 * (x / chord) - ...
        0.3516 * (x / chord).^2 + 0.2843 * (x / chord).^3 - 0.1036 * (x / chord).^4);
end