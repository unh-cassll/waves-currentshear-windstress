% Given wind speed, fetch, and wavenumber spectrum model parameter, apply
% the framework of Mueller & Veron [2009] to parse out:
% tangential viscous stress
% wave form drag
% stress associated with airflow separation from breaking wave crests
%
% Nathan Laxague 2021-2024
% in collaboration with David Ortiz-Suslow, Milan Curcic, and Jan-Victor Bjorkqvist
%
function out_struc = compute_stress_partition_MV09(air_side_friction_velocity_m_s,U10_m_s,air_density_kg_m3,wave_F_k_theta,wave_k_rad_m,wave_theta_rad,water_depth_m)

%% Parse input parameters

% Constants
g = 9.81;       % acceleration due to gravity in m/s^2
sigma = 0.072;  % surface tension in N/m
rho_w = 1030;   % water density in kg/m^3
rho_a = air_density_kg_m3;  % air density in kg/m^3
tau_total = rho_a*air_side_friction_velocity_m_s^2;

% Compute wave frequency and celerity given wavenumber
omega_rad_s = sqrt((g*wave_k_rad_m+sigma/rho_w*wave_k_rad_m.^3).*tanh(wave_k_rad_m*water_depth_m));
cp_m_s = omega_rad_s./wave_k_rad_m;

% Initial step: Plant [1982] parameterization
beta = 0.04*omega_rad_s.*(air_side_friction_velocity_m_s./cp_m_s).^2.*cos(wave_theta_rad);
integrand = rho_w*g*beta.*wave_F_k_theta./cp_m_s.*wave_k_rad_m.*cos(wave_theta_rad);
integrand(isnan(integrand)) = 0;
tau_wave_k = squeeze(trapz(wave_theta_rad,integrand,2));

[~,f2,f3,Usep,epsilon_b,Lambda_differential] = compute_airflow_separation_parameters(air_side_friction_velocity_m_s,U10_m_s,wave_k_rad_m,wave_theta_rad,wave_F_k_theta,beta);

tau_sep_int = get_stress_airflow_separation(f3,Usep,epsilon_b,Lambda_differential,air_density_kg_m3,wave_theta_rad,wave_k_rad_m,omega_rad_s);

tau_wave_k(isnan(tau_wave_k)) = 0;
tau_wave_k = abs(tau_wave_k);

% Separate viscous stress and form drag
tau_wave_int = trapz(wave_k_rad_m,tau_wave_k);
shorter_waves_drag = tau_wave_int - cumtrapz(wave_k_rad_m,tau_wave_k);
tau_visc_int = tau_total - tau_wave_int - tau_sep_int;

% Mueller & Veron [2009] approach
theta_inds = abs(wave_theta_rad) < pi;
Cb = 0.02 - 0.02*tanh(cp_m_s./(2*air_side_friction_velocity_m_s)-1.8*pi);
h_block = cos(wave_theta_rad).^1;
beta = 1/rho_a*Cb.*omega_rad_s.*cp_m_s.^-2.*h_block.*(tau_visc_int + tau_sep_int + shorter_waves_drag);
integrand = rho_w*g*beta(:,theta_inds).*wave_F_k_theta(:,theta_inds)./cp_m_s.*wave_k_rad_m.*cos(wave_theta_rad(theta_inds));
integrand(isnan(integrand)) = 0;
tau_wave_k = squeeze(trapz(wave_theta_rad(theta_inds),integrand,2));
tau_wave_k = abs(tau_wave_k);

% Feedback
tau_wave_int = trapz(wave_k_rad_m,f2.*tau_wave_k);
tau_sep_int = tau_total - tau_visc_int - tau_wave_int;
if tau_sep_int < 0
    tau_sep_int = 0;
end

% Output
out_struc.tau_total = tau_total;
out_struc.tau_visc = tau_visc_int;
out_struc.tau_wave = tau_wave_int;
out_struc.tau_sep = tau_sep_int;
out_struc.tau_wave_k = tau_wave_k;

end
