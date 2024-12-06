% Given air-side friction velocity, stress fraction, wind direction, and an
% observed subsurface current profile, this function stitches together the
% observed profile and the parameterized wind drift profile
%
% Nathan Laxague 2020-2024
%
function [z_current,U_current,dir_current_rel_wind,U_drift] = stitch_obs_to_wind_drift(air_side_friction_velocity_m_s_value,stress_fraction,wind_dir_rad,current_z_m,current_speed_m_s,current_dir_rad,wind_drift_param_flag)

if nargin<2
    stress_fraction = 1;
end

if nargin<3
    wind_dir_rad = 0;
end

if nargin<4
    current_z_m = -10;
end

if nargin<5
    current_speed_m_s = 0;
end

if nargin<6
    current_dir_rad = 0;
end

% Ensure that input profiles are in column form and go from shallow -> deep
current_z_m = current_z_m(:);
current_speed_m_s = current_speed_m_s(:);
current_dir_rad = current_dir_rad(:);
[current_z_m,order] = sort(current_z_m,'descend');
current_speed_m_s = current_speed_m_s(order);
current_dir_rad = current_dir_rad(order);

% Grab uppermost depth and current direction values
uppermost_obs_depth_m = max(current_z_m,[],'all','omitnan');

% Compute wind drift profile
[wind_drift_z,wind_drift_U] = compute_wind_drift_profile(air_side_friction_velocity_m_s_value,uppermost_obs_depth_m,stress_fraction);
[wind_drift_z,order] = sort(wind_drift_z,'descend');
wind_drift_U = wind_drift_U(order,:);
U_drift = max(wind_drift_U) - min(wind_drift_U);

% Decompose observed profile and parameterized wind drift into N & E
in_U_N = current_speed_m_s.*cos(current_dir_rad);
in_U_E = current_speed_m_s.*sin(current_dir_rad);
wind_drift_U_N = wind_drift_U*cos(wind_dir_rad);
wind_drift_U_E = wind_drift_U*sin(wind_dir_rad);

% Shift wind drift profile to be continuous with observation
wind_drift_U_N = wind_drift_U_N - wind_drift_U_N(end) + in_U_N(1);
wind_drift_U_E = wind_drift_U_E - wind_drift_U_E(end) + in_U_E(1);

% Concatenate profiles
num_winds = length(air_side_friction_velocity_m_s_value);
out_z = [wind_drift_z; current_z_m];
out_U_N = [wind_drift_U_N; repmat(in_U_N,1,num_winds)];
out_U_E = [wind_drift_U_E; repmat(in_U_E,1,num_winds)];

% Compute directional profiles in mag/dir format
z_current = out_z;
U_current = sqrt(out_U_N.^2+out_U_E.^2);
dir_current_rel_wind = mod(atan2(out_U_E,out_U_N)*180/pi - wind_dir_rad*180/pi,360)*pi/180;

if ~wind_drift_param_flag
    z_current = current_z_m;
    U_current = current_speed_m_s;
    dir_current_rel_wind = current_dir_rad;
    U_drift = 0;
end

