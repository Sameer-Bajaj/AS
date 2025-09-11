%{

Sameer Bajaj
Supersonic Airfoil Optimization Script

Computes drag impulse across leading and trailing edge using θ-β-M relation
for oblique shocks and prandtl-meyer function for expansion fans. Sources
are linked below each section. This program is only relevant for supersonic
flight, since waves are only generated in those regimes.

%}

clear
close all
clc

%% Import flight info
%{
Make sure your openrocket exports time, vertical velocity (ft/s), mach number,
reynolds number, temperature (F) & pressure (mbar). This is designed to
work w/ default imperical units
%}
% Insert your file name below
% Set the minimum mach to analyze for.
min_mach = 1.4;
flight_info = readmatrix("ares_flight.csv");
% deletes rows after recovery
flight_info = flight_info(flight_info(:,2) >= 0, :);
% removes non supersonic regions (lambda equation is best if mach > 1.2)
flight_info = flight_info(flight_info(:,3) >= min_mach, :);
% unit conversion

in2m = 0.0254;
ft2m = in2m*12;
N2lb = 1/4.4482216153;
mbar2Pa = 100;

flight_info(:, 2) = flight_info(:, 2)*ft2m;
flight_info(:, 4) = 5/9*(flight_info(:, 4) - 32) + 273.15;
flight_info(:, 5) = mbar2Pa * flight_info(:, 5);

time = flight_info(:, 1); %all rows, col 1
u_inf = flight_info(:,2); %flight velocity
mach = flight_info(:,3);%mach number
T = flight_info(:,4);
P_atm = flight_info(:,5);

%% Getting environmental factors
% PAGE 276: https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf
R = 8.314510; % univ gas const. in J/(mol-K)
a1 =  1.009950160D+04; a2 = -1.968275610D+02; a3 = 5.009155110D+00;
a4 = -5.761013730D-03; a5 = 1.066859930D-05; a6 = -7.940297970D-09;
a7 =  2.185231910D-12;
c_p = R*(a1*T.^-2 + a2*T.^-1 + a3 + a4*T + a5*T.^2 + a6*T.^3 + a7*T.^4);
gamma = c_p./(c_p - R);


%% Fin geometry
% fill this in based on your fins
c_tip = 8;
span = 8 * in2m;
t_half = .25 / 2; % half fin thickness
% min angle to analyze, the angle that has a wedge length equal to the tip
% chord
theta_min = rad2deg(atan(t_half/c_tip));
theta_max = 20; % maximum angle to analyze
thet = linspace(theta_min, theta_max, 15*theta_max);
% leading edge
alpha = deg2rad(thet);
% trailing edge
beta = alpha;

wedge_len = @(thet_in) t_half/tan(thet_in);

% preallocation
nt = numel(time); na = numel(alpha); nb = numel(beta);

D_1 = zeros(na, nt);                         % leading-edge drag [N] over time
D_2 = zeros(nb, nt);                         % trailing-edge drag [N] over time
ImpTot = nan(na, nb);                       % integrated impulse per (alpha,beta) pair [N*s]
Imp1 = D_1; Imp2 = D_2;
options = optimset("Display","off","FunValCheck","off");


for i = 1:na
    for j = 1:nb
        % set current values
        a = alpha(i); b = beta(j);
        len_lead = wedge_len(a)*in2m; % hypotenuse of leading wedge
        len_trail = wedge_len(b)*in2m; % "..." of trailing wedge
        % area of leading edge wedge
        A1 = span*len_lead;
        % area of trailing edge wedge
        A2 = span*len_trail;
        %% Solving θ-β-M relation:
        % https://web.archive.org/web/20121021100737/http://www.aerostudents.com/files/aerodynamicsC/obliqueShockWaves.pdf       
        x = linspace(deg2rad(1), pi/2, 500);
        for w = 1:nt
            M_1 = mach(w);
            Y = gamma(w);
            P_1 = P_atm(w);
            Fun = @(B) 2/tan(B)*(M_1^2*sin(B)^2-1)/(M_1^2*(Y+cos(2*B))+2)-tan(a);
            vec_Fun = 2./tan(x).*(M_1^2*sin(x).^2-1)./(M_1^2*(Y+cos(2*x))+2)-tan(a);
            % we need to get good bounds to find the first value of beta, which
            % will correspond to a weak shock. otherwise fzero could return a
            % higher root which we don't want.
            thet_2 = 0;
            for k = 1:length(x)-1
                % find the points between which we switch sign, and make sure
                % it's not an asymptote (delta < 100)
                if vec_Fun(k) * vec_Fun(k+1) <= 0 && abs(vec_Fun(k+1) - vec_Fun(k)) < 1e2
                    thet_2 = fzero(Fun, [x(k-1), x(k+1)], options);
                    % if at an asymptote, beta will return 0 from fzero. if
                    % this happens we wanna keep iterating
                    if thet_2 == 0
                        continue
                    else
                        break
                    end
                end
            end
            if thet_2 == 0 || not(isreal(thet_2))
                D_1(i, w) = NaN;
                continue
            end
            M_n_2 = sqrt(abs((2+(Y-1)*M_1^2*sin(thet_2)^2)/(2*Y*M_1^2*sin(thet_2)^2-(Y-1))));
            M_2 = M_n_2/sin(thet_2-a);
            P_2 = P_1*(1+2*Y/(Y+1)*(M_1^2*sin(thet_2)^2-1));
            % A1 * p1 is the force normal to the wedge, mult by sin(wedge 
            % angle) to get drag component. mult by 2 to account for lower
            % surface
            % D_1 contains info where row # corr. to ORK file step and col #
            % corr. to angle
            D_1(i, w) = 2*N2lb*P_2*A1*sin(a);
            if M_2 < 1
                D_1(i, w) = NaN;
                continue
            end
            %% Expansion fan!!
            % v_2 is the prandtl-meyer function at the second mach (after the
            % expansion fan)
            % https://www.grc.nasa.gov/www/k-12/airplane/pranmyer.html
            v = @(M_in) sqrt((Y+1)/(Y-1))*atan(sqrt((Y-1)/(Y+1)*(M_in^2-1))) ...
                - atan(sqrt(M_in^2-1));
            v_2 = v(M_2);
            fun_2 = @(M_in) v(M_in) - v_2 - a;
            M_3 = fzero(fun_2, M_2, options); % find the new Mach number after the expansion fan at leading edge
            if M_3 <= 1 || isnan(M_3)
                D_1(i,w) = NaN;
                continue
            end
            v_3 = v(M_3);
            P_3 = P_2*((1+(Y-1)/2*M_2^2)/(1+(Y-1)/2*M_3^2))^(Y/(Y-1));
            fun_3 = @(M_in) v(M_in) - v_3 - b;
            %% Given this P_3, calculate the expansion fan drag for the trailing edge.
            if not(isreal(v_3)) || v_3 == 0
                D_1(i,w) = NaN;
                continue
            end
            M_4 = fzero(fun_3, M_3, options);
            if M_4 < 1
                D_2(j, w) = NaN;
                continue
            end
            P_4 = P_3*((1+(Y-1)/2*M_3^2)/(1+(Y-1)/2*M_4^2))^(Y/(Y-1));
            D_2(j, w) = 2*N2lb*P_4*A2*sin(b); % calculate drag for the trailing edge wedge
        end
        Imp1(i) = trapz(time, D_1(i, :));
        Imp2(j) = trapz(time, D_2(j, :));
        ImpTot(i, j) = Imp1(i) + Imp2(j); % N*s
    end
end

flight_info = flight_info(flight_info(:,2) >= 0, :);
[a_safe, ~] = find(isnan(D_1), 1, 'first');
[b_safe, ~] = find(isnan(D_2), 1, 'first');
% leading edge drag vs angle
figure
plot(thet, Imp1, 'LineWidth', 1)
xlim([theta_min thet(a_safe)])
xlabel('Leading Edge Wedge Angle (deg)')
ylabel('Drag Impulse [N\bullets]')
title('Leading Edge Drag Impulse vs \alpha')
grid on
% TE drag v angle
figure
plot(thet, Imp2, 'LineWidth', 1)
xlim([theta_min thet(b_safe)])
xlabel('Trailing Edge Wedge Angle (deg)')
ylabel('Drag Impulse [N\bullets]')
title('Trailing Edge Drag Impulse vs \beta')
grid on
% LE drag v length
wedge_len_vec = t_half./tan(deg2rad(thet));
figure
plot(wedge_len_vec, Imp1, 'LineWidth', 1)
xlim([wedge_len_vec(a_safe) wedge_len(deg2rad(theta_min))])
xlabel('Leading Edge Wedge Len (in)')
ylabel('Drag Impulse [N\bullets]')
title('Drag Impulse vs Leading Edge Length')
grid on
% TE drag v length
figure
plot(wedge_len_vec, Imp2, 'LineWidth', 1)
xlim([wedge_len_vec(b_safe) wedge_len(deg2rad(theta_min))])
xlabel('Trailing Edge Wedge Len (in)')
ylabel('Drag Impulse [N\bullets]')
title('Drag Impulse vs Trailing Edge Length')
grid on

[BetaGrid, AlphaGrid] = meshgrid(thet, thet);
[LA, LB] = ndgrid(wedge_len_vec, wedge_len_vec);  % size na×nb
mask = (LA + LB) > c_tip;
ImpTot(mask) = NaN;
figure
surf(AlphaGrid, BetaGrid, ImpTot, 'EdgeColor', 'none')
xlabel('\alpha [deg]')
ylabel('\beta [deg]')
zlabel('Drag Impulse (N\bullets)')
title('Wedge Drag Impulse vs Leading / Trailing Wedge Angle')
colorbar
view(135, 30)

[~, lin_idx] = min(ImpTot, [], 'all');
alpha_idx = mod(lin_idx, size(ImpTot, 1));
beta_idx = round((lin_idx-alpha_idx)/size(ImpTot, 2));
fprintf('Optimal LE angle: %.4f deg\nOptimal TE angle: %.4f deg\n', ...
    thet(alpha_idx), thet(beta_idx));