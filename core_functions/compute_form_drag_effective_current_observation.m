% Given wind stress, the observed current profile, and the wave spectrum,
% computes the normalized form drag, effective current speed, and effective
% current depth
%
% Nathan Laxague 2024
% in collaboration with David Ortiz-Suslow, Milan Curcic, and Jan-Victor Bjorkqvist
%
function out_struc = compute_form_drag_effective_current_observation(obs_wind,obs_current_profile,obs_wave_spectrum,water_depth_m,k_max,wind_drift_param_flag)

% Form of input parameters

% form of 'obs_wind':
% * air density in kg/m^3 -> obs_wind.air_density_kg_m3
% ** single value
% * air-side friction velocity in m/s -> obs_wind.air_side_friction_velocity_m_s
% ** single value
% * wind speed in m/s -> obs_wind.wind_speed_m_s
% ** single value
% * wind direction (coming-from) in rad -> obs_wind.wind_direction_rad
% ** single value
% * stress angle in rad -> obs_wind.stress_angle_rad
% ** single value
% * wind measurement height in m -> obs_wind.wind_meas_height_m
% ** single value

% form of 'obs_current_profile':
% * depth [-depth 0] in m -> obs_current_profile.z_m
% ** Nx1 array
% * speed in m/s -> obs_current_profile.current_speed_m_s
% ** Nx1 array
% * direction going-to in rad -> obs_current_profile.current_dir_rad
% ** Nx1 array

% form of 'obs_wave_spectrum':
% * frequency in Hz -> obs_wave_spectrum.f_Hz
% ** Mx1 array
% * direction going-to in rad -> obs_wave_spectrum.theta_rad
% ** 1xL array
% * spectral energy density in m^2/Hz/rad -> obs_wave_spectrum.F_f_theta
% ** MxL array

% form of 'water_depth_m':
% * water depth in m
% ** single value

% Wind
air_density_kg_m3 = obs_wind.air_density_kg_m3;
air_side_friction_velocity_m_s = obs_wind.air_side_friction_velocity_m_s;
stress_angle_rad = obs_wind.stress_angle_rad;
wind_dir_rad = pi/180*mod((obs_wind.wind_direction_rad + stress_angle_rad)*180/pi+180,360);
wind_speed_m_s_value = obs_wind.wind_speed_m_s;
wind_meas_height_m = obs_wind.wind_meas_height_m;
U10_m_s = air_side_friction_velocity_m_s/0.4*log(10/wind_meas_height_m) + wind_speed_m_s_value;

% Current profile
current_z_m = obs_current_profile.z_m;
current_speed_m_s = obs_current_profile.current_speed_m_s;
current_dir_rad = obs_current_profile.current_dir_rad;

% Wave spectrum
wave_f_Hz = obs_wave_spectrum.f_Hz;
wave_theta_rad = obs_wave_spectrum.theta_rad;
wave_F_f_theta = obs_wave_spectrum.F_f_theta;

% Gather variable dimensions
num_wave_directions = length(wave_theta_rad);

% Transform to downwind/crosswind reference

% Wave spectrum
wave_F_f_theta_big = [wave_F_f_theta(:,1:end-1) wave_F_f_theta(:,1:end-1) wave_F_f_theta(:,1:end-1)];
wave_theta_rad_big = [wave_theta_rad(1:end-1)-2*pi wave_theta_rad(1:end-1) wave_theta_rad(1:end-1)+2*pi];
wave_theta_rad_big_downwind = wave_theta_rad_big - mod((wind_dir_rad*180/pi),360)*pi/180;
ind_start = find(wave_theta_rad_big_downwind>-pi,1,'first');
inds_to_consider = ind_start:ind_start-1+num_wave_directions;
wave_F_f_theta_rel_wind = wave_F_f_theta_big(:,inds_to_consider);
wave_theta_rad_rel_wind = wave_theta_rad_big((1:length(inds_to_consider))+floor(length(inds_to_consider)/2));

% Compute stress (and constituent components)

% viscous stress fraction starting guess
stress_fraction = 1;

% 1. given stress fraction and wind stress, parameterize wind drift
[current_z_m_with_wind_drift,current_speed_m_s_with_wind_drift,current_dir_deg_rel_wind,U_drift] = stitch_obs_to_wind_drift(air_side_friction_velocity_m_s,stress_fraction,wind_dir_rad,current_z_m,current_speed_m_s,current_dir_rad,wind_drift_param_flag);

% 2. compute wavenumber-directional spectrum given frequency spectrum and
% current profile
[wave_F_k_theta_rel_wind,wave_k_rad_m,U_current_k,D_current_k] = directional_Doppler_shift_spectrum(current_speed_m_s_with_wind_drift,current_dir_deg_rel_wind,current_z_m_with_wind_drift,water_depth_m,wave_f_Hz,wave_F_f_theta_rel_wind,wave_theta_rad_rel_wind,k_max);

% 3. send wind stress and wavenumber-directional spectrum to MV2009 stress
% partition function
stress_partition_struc = compute_stress_partition_MV09(air_side_friction_velocity_m_s,U10_m_s,air_density_kg_m3,wave_F_k_theta_rel_wind,wave_k_rad_m,wave_theta_rad_rel_wind,water_depth_m);

% Compute U_w and U_v

% U_w - weighted by form stress density
weighting_function = stress_partition_struc.tau_wave_k;
U_w = trapz(wave_k_rad_m,weighting_function.*U_current_k)./trapz(wave_k_rad_m,weighting_function);

% U_v - surface current
U_v = current_speed_m_s_with_wind_drift(1);

% Output
out_struc.U10_m_s = U10_m_s;
out_struc.k_rad_m = wave_k_rad_m;
out_struc.wave_F_k = trapz(wave_theta_rad_rel_wind,wave_k_rad_m.*wave_F_k_theta_rel_wind,2);
out_struc.tau_total = stress_partition_struc.tau_total;
out_struc.tau_w_k = stress_partition_struc.tau_wave_k;
out_struc.tau_w = stress_partition_struc.tau_wave;
out_struc.tau_v = stress_partition_struc.tau_visc;
out_struc.tau_s = stress_partition_struc.tau_sep;
out_struc.U_k = U_current_k;
out_struc.D_k = D_current_k;
out_struc.U_w = U_w;
out_struc.U_v = U_v;
out_struc.U_drift = U_drift;
out_struc.current_z_m = current_z_m_with_wind_drift;
out_struc.current_U_m_s = current_speed_m_s_with_wind_drift;
out_struc.current_D_deg = current_dir_deg_rel_wind;
