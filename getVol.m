function vol = getVol(len, OD, thickness)
    C = 0; % for von kaarman, C = 0, for LV-haack, C = 1/3
    t1 = 0; t2 = pi;
    
    y_in = @(t) ~(thickness == 0) * (OD-2*thickness)/2/sqrt(pi)*sqrt(t-0.5*sin(2*t)+C^3*sin(t).^3);
    y_out = @(t) (OD)/2/sqrt(pi)*sqrt(t-0.5*sin(2*t)+C^3*sin(t).^3);
    dx_dt_in = @(t) ~(thickness == 0) * (len-thickness)/2*(sin(t));
    dx_dt_out = @(t) (len)/2*sin(t);
    vol_in = ~(thickness == 0) * pi * integral(@(t) y_in(t).^2 .* dx_dt_in(t), t1, t2);
    vol_out = pi * integral(@(t) y_out(t).^2 .* dx_dt_out(t), t1, t2);
    vol = vol_out - vol_in; % in^3
end