%% Airfoil Optimization Script

clear
close all
clc

%% Import flight info
%{
Make sure your openrocket exports altitude, vertical velocity, mach number,
reynolds number, temperature & pressure
%}
% Insert your file name below
flight_info = readmatrix("ares_flight.csv");
% deletes rows after recovery
flight_info = flight_info(flight_info(:,2) >= 0, :);
% removes non supersonic regions (lambda equation is best if mach > 1.2)
flight_info = flight_info(flight_info(:,3) >= 1.2, :);
altitude = flight_info(:, 1); %all rows, col 1
u_inf = flight_info(:,2); %flight velocity
mach = flight_info(:,3);%mach number
Re_L = flight_info(:,4); %reynolds number
T_atm = flight_info(:,5);
P_atm = flight_info(:,6);

in2m = 0.0254;
ft2m = in2m*12;
N2lb = 1/4.4482216153;

%% Fin geometry
c_root = 18 * in2m;
c_tip = 8 * in2m;
span = 8 * in2m;
t_half = .25 * in2m / 2; % half fin thickness
% checks if wedge is gonna fit on the tip chord
theta_min = ceil(rad2deg(atan(2*t_half/c_tip)));
theta_max = 20; % maximum angle to analyze
theta = linspace(theta_min, theta_max, 15*(theta_max-theta_min));
rad = deg2rad(theta);
lambda = sqrt(mach.^2-1);
rho = P_atm./(287.058 * T_atm);
q_inf = 1/2*u_inf.^2.*rho;
A = sqrt((c_root-c_tip)^2+span^2) * t_half./sin(rad);
% initialize drag matrix
D = zeros(length(altitude), length(rad));

for i = 1:length(altitude)
    for j = 1:length(rad)
        len_wedge = t_half / tan(rad(j)); % horizontal width of wedge
        hyp_wedge = t_half / sin(rad(j)); % total length of wedge piece
        % area of leading edge wedge
        A1 = span*hyp_wedge;
        % area of trailing edge wedge
        A2 = span*hyp_wedge;
        % coefficient of pressure at specific mach and angle
        c_p = 2*rad(j)/lambda(i);
        % pressures in 1. fore and 2. aft wedges at certain altitude
        % https://www.sciencedirect.com/topics/engineering/critical-pressure-coefficient

        p1 = P_atm(i) + c_p*q_inf(i);
        % angle is negative for second part, so c_p (2) becomes negative
        p2 = P_atm(i) - c_p*q_inf(i);

        % A1 * p1 is the force normal to the wedge, mult by sin(wedge 
        % angle) to get drag component

        D(i, j) = 2*(A1*p1 - A2*p2)*sin(rad(j))*N2lb;
        % mult by 2 to account for bottom surface
    end
end

[ThetaGrid, AltGrid] = meshgrid(theta, altitude);

figure
surf(AltGrid, ThetaGrid, D, 'EdgeColor', 'none')
xlabel('Altitude [m]')
ylabel('\theta [deg]')
zlabel('Drag [lbf]')
title('Drag vs Altitude and Wedge Angle')
colorbar
view(135, 30) % nice angle



