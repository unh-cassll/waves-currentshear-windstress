% Computes separation stress
% following Mueller & Veron [2009]
%
% Code by N. Laxague 2024
%
function out_tau_sep_int = get_stress_airflow_separation(f3,Usep,epsilon_b,Lambda_differential,rho_a,wave_theta_rad,wave_k_rad_m,omega_rad_s)

gamma = 0.45;

cp_m_s = omega_rad_s./wave_k_rad_m;

Us = Usep.*cos(wave_theta_rad) - cp_m_s;

k_inds = wave_k_rad_m < 20*pi;
theta_inds = abs(wave_theta_rad)<pi;

out_tau_sep_int = rho_a*epsilon_b*gamma*trapz(wave_k_rad_m(k_inds),f3(k_inds).*trapz(wave_theta_rad(theta_inds),Us(k_inds,theta_inds).^2.*cos(wave_theta_rad(theta_inds)).*Lambda_differential(k_inds,theta_inds),2));

end