% Computes airflow separation parameters
% following Mueller & Veron [2009]
%
% Code by N. Laxague 2024
%
function [f1,f2,f3,Usep,epsilon_b,Lambda_differential] = compute_airflow_separation_parameters(air_side_friction_velocity_m_s,U10_m_s,wave_k_rad_m,wave_theta_rad,wave_F_k_theta,beta)

% Constants
g = 9.81;       % acceleration due to gravity in m/s^2
sigma = 0.072;  % surface tension in N/m
rho_w = 1030;   % water density in kg/m^3

omega_rad_s = sqrt(g*wave_k_rad_m+sigma/rho_w*wave_k_rad_m.^3);

B_rad = squeeze(trapz(wave_theta_rad,wave_k_rad_m.^4.*wave_F_k_theta,2));

epsilon_b = get_breaking_wave_slope(wave_k_rad_m,B_rad);
breaking_b = get_normalized_wave_breaking_dissipation_rate(epsilon_b);

L = (epsilon_b.^(1/2)+1/4);

L_mat = L*ones(length(wave_k_rad_m),length(wave_theta_rad)).*cos(wave_theta_rad);
breaking_b_mat = breaking_b*ones(length(wave_k_rad_m),length(wave_theta_rad));

L_tilde = beta.*wave_k_rad_m.^4.*wave_F_k_theta.*(omega_rad_s.*breaking_b_mat).^-1;
Pbr = 2*pi./wave_k_rad_m.*L_tilde;

Lambda_differential = L_tilde./wave_k_rad_m;

Q = 1 - squeeze(trapz(wave_theta_rad,L_mat.*Pbr,2));

A_tilde_integral = squeeze(trapz(wave_theta_rad,L_mat.*Pbr,2));

product_mat = cumprod(Q);

A_tilde = product_mat.*A_tilde_integral;

f1 = 1 - trapz(wave_k_rad_m,A_tilde);
f2 = 1 - cumtrapz(wave_k_rad_m,A_tilde);
f3 = A_tilde.*(squeeze(trapz(wave_theta_rad,Pbr,2))).^-1;

f1(isnan(f1)) = 0;
f2(isnan(f2)) = 0;
f3(isnan(f3)) = 0;

f2(f2<0) = 0;
f3(f3<0) = 0;

Usep = retrieve_Useparation_breaking(air_side_friction_velocity_m_s,U10_m_s,epsilon_b,wave_k_rad_m);

    function epsilon_b = get_breaking_wave_slope(wave_k_rad_m,B_rad)

        f_Hz = sqrt(g*wave_k_rad_m+sigma/rho_w*wave_k_rad_m.^3)/(2*pi);
        cg_m_s = 0.5*2*pi*f_Hz./wave_k_rad_m;

        S = wave_k_rad_m.^-2.*B_rad;

        Sf = 2*pi*wave_k_rad_m./cg_m_s.*S;

        ind_peak = find(Sf==max(Sf,[],'all','omitnan'),1,'first');
        ind_low = find(f_Hz>0.5*f_Hz(ind_peak),1,'first');
        ind_high = find(f_Hz>1.5*f_Hz(ind_peak),1,'first');
        inds = ind_low:ind_high;

        epsilon_b = 2*sqrt(trapz(wave_k_rad_m(inds),S(inds)));

    end

    function breaking_b = get_normalized_wave_breaking_dissipation_rate(epsilon_b)

        Upsilon = 0.01;
        chi = 0.25;

        inds_spill = epsilon_b <= 0.2;
        inds_plunge = epsilon_b > 0.2;

        breaking_b = NaN*epsilon_b;
        breaking_b(inds_spill) = Upsilon*epsilon_b(inds_spill).^(1/2);
        breaking_b(inds_plunge) = chi*epsilon_b(inds_plunge).^(5/2);

    end

    function [out_Uwind,out_zwind] = get_wind_profile(air_side_friction_velocity_m_s,U10_m_s)

        kappa = 0.4;

        out_zwind = (0:0.001:100)';
        out_Uwind = air_side_friction_velocity_m_s/kappa*log(out_zwind/10) + U10_m_s;
        out_Uwind(out_Uwind<0) = 0;

    end

    function Usep = retrieve_Useparation_breaking(air_side_friction_velocity_m_s,U10_m_s,epsilon_b,wave_k_rad_m)

        [out_Uwind,out_zwind] = get_wind_profile(air_side_friction_velocity_m_s,U10_m_s);

        zsep = epsilon_b(:)./wave_k_rad_m(:);
        Usep = NaN*zsep;

        for i = 1:length(wave_k_rad_m)

            zbit = zsep(i);
            zdiff = abs(zbit-out_zwind);
            ind = find(zdiff==min(zdiff,[],'all','omitnan'),1,'first');
            Usep(i) = out_Uwind(ind);

        end

    end

end