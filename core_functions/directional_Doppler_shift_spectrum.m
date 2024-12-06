% Given current profile and frequency-directional spectrum, computs
% wavenumber-directional spectrum via Doppler shift
%
% Nathan Laxague 2020-2024
%
function [wave_F_k_theta_rel_wind,k_rad_m,U_current_k,D_current_k] = directional_Doppler_shift_spectrum(current_speed_m_s_with_wind_drift,current_dir_deg_rel_wind,current_z_m_with_wind_drift,water_depth_m,wave_f_Hz,wave_F_f_theta_rel_wind,wave_theta_rad_rel_wind,k_max)

% Constants
g = 9.81;       % acceleration due to gravity in m/s^2
sigma = 0.072;  % surface tension in N/m
rho_w = 1030;   % water density in kg/m^3

% Compute U(k) from U(z) - Equation 18 from Kirby & Chen [1989], JGR
% Input arrays of size NxM (depths by frequencies)
% Output array U(k) of size Mx1 (number of frequencies/wavenumbers)
k_reference = logspace(-4,log10(5),256)';
wave_theta_rad_rel_wind = wave_theta_rad_rel_wind(:)';
wave_f_Hz = wave_f_Hz(:);
[s1,~] = size(wave_F_f_theta_rel_wind);
if length(wave_theta_rad_rel_wind) == s1
    wave_F_f_theta_rel_wind = wave_F_f_theta_rel_wind';
end

num_current_depths = length(current_z_m_with_wind_drift);

U_block = repmat(current_speed_m_s_with_wind_drift(:),[1 length(k_reference)]);
D_block = repmat(current_dir_deg_rel_wind(:),[1 length(k_reference)]);
z_block = repmat(current_z_m_with_wind_drift(:),[1 length(k_reference)]);

% Compute U(k) given U(z)
k_block = repmat(k_reference(:)',[num_current_depths 1]);
front_part = 2*k_block./(sinh(2*k_block*water_depth_m));

integrand_downwind = U_block.*cosd(D_block).*cosh(2*k_block.*(water_depth_m+z_block));
U_current_k_downwind = -1*squeeze(trapz(current_z_m_with_wind_drift,front_part.*integrand_downwind))';
U_current_k_downwind(isnan(U_current_k_downwind)) = 0;

integrand_crosswind = U_block.*sind(D_block).*cosh(2*k_block.*(water_depth_m+z_block));
U_current_k_crosswind = -1*squeeze(trapz(current_z_m_with_wind_drift,front_part.*integrand_crosswind))';
U_current_k_crosswind(isnan(U_current_k_crosswind)) = 0;

U_current_k = sqrt(U_current_k_downwind.^2+U_current_k_crosswind.^2);
D_current_k = atan2(U_current_k_crosswind,U_current_k_downwind);

% Compute reference wave frequency array from a test wavenumber and a
% known current velocity
f_mat = sqrt((g*k_reference+sigma/rho_w*k_reference.^3).*tanh(k_reference*water_depth_m)+k_reference.*U_current_k.*cos(D_current_k-wave_theta_rad_rel_wind))/(2*pi);

% Find matching wavenumber (numerical inversion of LDR)
k_rad_m = NaN*wave_F_f_theta_rel_wind;
for i = 1:length(wave_f_Hz)
    for j = 1:length(wave_theta_rad_rel_wind)
        f_diff = abs(f_mat(:,j)-wave_f_Hz(i));
        ind = find(f_diff==min(f_diff),1,'first');
        k_rad_m(i,j) = k_reference(ind);
    end
end

% Convert from F(f,theta) to F(k,theta), conserving energy
Cg = 2*pi*wave_f_Hz./k_rad_m/2.*(1+2*k_rad_m*water_depth_m./sinh(2*k_rad_m*water_depth_m));
wave_F_k_theta_rel_wind_converted = wave_F_f_theta_rel_wind.*Cg./(2*pi*k_rad_m);

% Interpolates onto uniform wavenumber grid
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
k_vec = linspace(max(k_rad_m(1,:)),min(k_rad_m(end,:)),128)';
F = scatteredInterpolant(reshape(k_rad_m,[],1),reshape(repmat(wave_theta_rad_rel_wind,length(wave_f_Hz),1),[],1),reshape(wave_F_k_theta_rel_wind_converted,[],1));
wave_F_k_theta_rel_wind_interp_vec = F(reshape(repmat(k_vec,1,length(wave_theta_rad_rel_wind)),[],1),reshape(repmat(wave_theta_rad_rel_wind,length(k_vec),1),[],1));
wave_F_k_theta_rel_wind_interp_mat = reshape(wave_F_k_theta_rel_wind_interp_vec,length(k_vec),length(wave_theta_rad_rel_wind));

% Attach k^-3 tail (including to directional spreading function)
in_struc.k_rad_m = k_vec;
in_struc.theta_rad = wave_theta_rad_rel_wind;
in_struc.F_k_theta = wave_F_k_theta_rel_wind_interp_mat;
tail_struc = pin_the_tail_on_the_spectrum(in_struc,k_max);
k_rad_m = tail_struc.k_rad_m;
wave_F_k_theta_rel_wind = tail_struc.F_k_theta_new;

U_current_k = interp1(k_reference,U_current_k,k_vec);
D_current_k = interp1(k_reference,D_current_k,k_vec);