%
function [out_z,out_U] = compute_wind_drift_profile(u_star_a,max_depth,stress_fraction)

% Ensure that 'ustar_a' = row, 'stress_fraction' = column, 'max_depth' < 0
u_star_a = u_star_a(:)';
stress_fraction = stress_fraction(:)';
max_depth = -1*abs(max_depth);

% Compute vertical coordinate vector
dz = 0.01;
zeta = (max_depth+dz:dz:-dz)';

% Constants
rho_a = 1.226;      % Density of air, kg/m^3
rho_w = 1030;       % Density of water, kg/m^3
visc_water = 1e-6;	% Kinematic viscosity of water, m^2/s

% Compute air and water-side stresses, 
tau_a = rho_a*u_star_a.^2;
tau_w = stress_fraction.*tau_a;
ustar_w = sqrt(tau_w/rho_w);

% Convert to wall coordinates
zeta_plus = abs(zeta)*ustar_w*visc_water^-1;
u_plus = zeta_plus;
inds_vsl = zeta_plus <= 5;
inds5 = zeta_plus > 5;
inds30 = zeta_plus > 30;

% Enforce piecewise profile described by Spalding (1961)
u_plus(inds_vsl) = zeta_plus(inds_vsl);
u_plus(inds5) = 2/0.4*log(zeta_plus(inds5)) - 3.05;
u_plus(inds30) = 1/0.4*log(zeta_plus(inds30)) + 5.5;

% Convert back
U_zeta_wall = u_plus.*repmat(ustar_w,length(zeta),1);
U_max = repmat(nanmax(U_zeta_wall),length(zeta),1);
U_zeta = U_max - U_zeta_wall;

out_z = zeta;
out_U = U_zeta;
