
% This script was run to ingest OOI/NDBC data collected from the
% Coastal Endurance Oregon Shelf location as part of the manuscript
% "Accounting for Waves and Current Shear in Ocean Wind Stress Parameterization"
% 
% D. G. Ortiz-Suslow, N. J. M. Laxague, M. Curcic, & J.-V. I. Bjorkqvist, 2024
%

% Outputs from this script are contained within the file
% '2024_EFM_manuscript/_data/OOI_CE_data_compilation.mat'

% Code written by Nathan Laxague (2024)

close all;clear;clc

load('../ooi_data/Coastal_Endurance/oregon_shelf/CE_oregon_shelf_bulk_met.mat')
load('../ooi_data/Coastal_Endurance/oregon_shelf/CE_oregon_shelf_direct_covariance_fluxes.mat')
load('../ooi_data/Coastal_Endurance/oregon_shelf/CE_oregon_shelf_single_point_velocity.mat')
load('../ooi_data/Coastal_Endurance/oregon_shelf/CE_oregon_shelf_ADCP.mat')
NDBC_struc = load('../ndbc_data/CE_oregon_shelf_NDBC46097.mat');

% three hour moving average
crosswind_covariance_m2_s2 = smoothdata(crosswind_covariance_m2_s2,'movmean',3,'omitnan');
downwind_covariance_m2_s2 = smoothdata(downwind_covariance_m2_s2,'movmean',3,'omitnan');
single_point_vN = smoothdata(single_point_vN,'movmean',13,'omitnan');
single_point_vE = smoothdata(single_point_vE,'movmean',13,'omitnan');

ADCP_z = -55:1:-10;
ADCP_DN = datenum(DTime_600kHz_ADCP(1)):15*60/86400:datenum(DTime_600kHz_ADCP(end));
ADCP_DTime = datetime(datevec(ADCP_DN),'TimeZone','UTC');
ADCP_E = NaN*ones(length(ADCP_z),length(ADCP_DN));
ADCP_N = ADCP_E;
ADCP_U = ADCP_E;
for n = 1:length(ADCP_z)

    inds = find(ADCP_bin_depth_m>ADCP_z(n)-0.5 & ADCP_bin_depth_m<ADCP_z(n)+0.5);
    DTime_particular = DTime_600kHz_ADCP(inds);
    ADCP_e = ADCP_eastward_sea_water_velocity(inds);
    ADCP_n = ADCP_northward_sea_water_velocity(inds);
    ADCP_u = ADCP_upward_sea_water_velocity(inds);

    ADCP_e = smoothdata(ADCP_e,'movmean',13,'omitnan');
    ADCP_n = smoothdata(ADCP_n,'movmean',13,'omitnan');
    ADCP_u = smoothdata(ADCP_u,'movmean',13,'omitnan');

    ADCP_E(n,:) = interp1(DTime_particular,ADCP_e,ADCP_DTime,'linear');
    ADCP_N(n,:) = interp1(DTime_particular,ADCP_n,ADCP_DTime,'linear');
    ADCP_U(n,:) = interp1(DTime_particular,ADCP_u,ADCP_DTime,'linear');

end

ADCP_inds_remove = ADCP_z <=-15 & ADCP_z >= -18;
ADCP_E = interp1(ADCP_z(~ADCP_inds_remove),ADCP_E(~ADCP_inds_remove,:),ADCP_z,'linear');
ADCP_N = interp1(ADCP_z(~ADCP_inds_remove),ADCP_N(~ADCP_inds_remove,:),ADCP_z,'linear');
ADCP_U = interp1(ADCP_z(~ADCP_inds_remove),ADCP_U(~ADCP_inds_remove,:),ADCP_z,'linear');

ADCP_E = smoothdata2(ADCP_E,'movmean',{3,1});
ADCP_N = smoothdata2(ADCP_N,'movmean',{3,1});
ADCP_U = smoothdata2(ADCP_U,'movmean',{3,1});

water_depth_m = 80;

k_max = 50;
f_cut = 0.35;

shallow_depth = -0.1;
filt_num = 11;

in_folder = '../ooi_data/Coastal_Endurance/oregon_shelf/wave_spectra/';

yearvec = 2017:2023;

out_path = '../output_data/';

for year_ind = 1:length(yearvec)


    file_info = dir([in_folder '*' num2str(yearvec(year_ind)) '*.nc']);

    for file_num = 1:length(file_info)

        nc_name = [in_folder file_info(file_num).name];

        DV = ncread(nc_name,'DateVector');
        DTime_vec = datetime(DV,'TimeZone','UTC');

        theta_rad_vec = pi/180*double(ncread(nc_name,'theta_deg'));

        f_Hz_vec = double(ncread(nc_name,'f_Hz_DIRSPEC'));
        F_f_theta_block = 10.^double(ncread(nc_name,'log10_F_f_theta_DIRSPEC'));

        F_f_theta_block = medfilt3(F_f_theta_block,[1 1 3]);

        num_spectra = size(F_f_theta_block,3);

        month_struc = struct();

        parfor spec_num = 1:num_spectra

            DTime = DTime_vec(spec_num);

            F_f_theta = squeeze(F_f_theta_block(:,:,spec_num));

            if sum(F_f_theta,'all') > 0

                DN_diff_wind = abs(datenum(DTime)-datenum(NDBC_struc.DTime));
                DN_diff_fluxes = abs(datenum(DTime)-datenum(DTime_fluxes));
                DN_diff_bulk = abs(datenum(DTime)-datenum(DTime_bulk));
                DN_diff_single_point = abs(datenum(DTime)-datenum(DTime_single_point_velocity));
                DN_diff_ADCP = abs(datenum(DTime)-datenum(ADCP_DTime));
                DN_diff_HFR_2km = abs(datenum(DTime)-datenum(HFR_2km_DTime));
                DN_diff_HFR_6km = abs(datenum(DTime)-datenum(HFR_6km_DTime));

                inds_wind = find(DN_diff_wind<1/24);
                inds_fluxes = find(DN_diff_fluxes<1/24);
                inds_bulk = find(DN_diff_bulk<1/24);
                inds_single_point = find(DN_diff_single_point<1/24);
                inds_ADCP = find(DN_diff_ADCP<1/24);
                inds_HFR_2km = find(DN_diff_HFR_2km<1/24);
                inds_HFR_6km = find(DN_diff_HFR_6km<1/24);

                wind_flag = false;
                single_point_flag = false;
                ADCP_flag = false;
                HFR_2km_flag = false;
                HFR_6km_flag = false;

                if ~isempty(inds_wind)
                    wind_dir_deg_coming_from = mod(180/pi*meanang(NDBC_struc.wind_dir_coming_from_deg(inds_wind)*pi/180)+360,360);
                    windspeed_m_s = median(NDBC_struc.wind_speed_m_s(inds_wind),'omitnan');
                    air_density_kg_m3_value = median(air_density_kg_m3(inds_bulk),'omitnan');
                    downwind_stress_N_m2 = -air_density_kg_m3_value*median(downwind_covariance_m2_s2(inds_fluxes),'omitnan');
                    crosswind_stress_N_m2 = -air_density_kg_m3_value*median(crosswind_covariance_m2_s2(inds_fluxes),'omitnan');
                    ustar_m_s = (downwind_stress_N_m2^2+crosswind_stress_N_m2^2)^(1/4)*air_density_kg_m3_value^(-1/2);
                    stress_angle_rad = atan2(crosswind_stress_N_m2,downwind_stress_N_m2);
                    if ~isnan(windspeed_m_s) && ~isnan(stress_angle_rad)
                        wind_flag = true;
                    end
                end

                if ~isempty(inds_single_point)
                    single_point_depth_m = median(single_point_velocity_depth_m(inds_single_point),'omitnan');
                    single_point_E = median(single_point_vE(inds_single_point),'omitnan');
                    single_point_N = median(single_point_vN(inds_single_point),'omitnan');
                    single_point_dir_deg_going_to = mod(atan2(single_point_E,single_point_N)*180/pi+360,360);
                    single_point_speed_m_s = sqrt(single_point_E^2+single_point_N^2);
                    if ~isnan(single_point_speed_m_s)
                        single_point_flag = true;
                    end
                    z_m = fliplr(-1*logspace(log10(single_point_depth_m),log10(30),1000))';
                end

                if ~isempty(inds_ADCP)
                    ADCP_E_median = median(ADCP_E(:,inds_ADCP),2,'omitnan');
                    ADCP_N_median = median(ADCP_N(:,inds_ADCP),2,'omitnan');
                    ADCP_speed_m_s = sqrt(ADCP_E_median.^2+ADCP_N_median.^2);
                    ADCP_dir_deg_going_to = mod(atan2(ADCP_E_median,ADCP_N_median)*180/pi+360,360);
                    if ~isnan(mean(ADCP_speed_m_s,'omitnan'))
                        ADCP_flag = true;
                    end
                end

                if ~isempty(inds_HFR_2km)
                    HFR_2km_E = median(HFR_2km_eastward_sea_water_velocity_m_s(inds_HFR_2km),'omitnan');
                    HFR_2km_N = median(HFR_2km_northward_sea_water_velocity_m_s(inds_HFR_2km),'omitnan');
                    HFR_2km_speed_m_s = sqrt(HFR_2km_E^2+HFR_2km_N^2);
                    HFR_2km_dir_deg_going_to = mod(atan2(HFR_2km_E,HFR_2km_N)*180/pi+360,360);
                    if ~isnan(HFR_2km_speed_m_s)
                        HFR_2km_flag = true;
                    end
                end

                if ~isempty(inds_HFR_6km)
                    HFR_6km_E = median(HFR_6km_eastward_sea_water_velocity_m_s(inds_HFR_2km),'omitnan');
                    HFR_6km_N = median(HFR_6km_northward_sea_water_velocity_m_s(inds_HFR_2km),'omitnan');
                    HFR_6km_speed_m_s = sqrt(HFR_6km_E^2+HFR_6km_N^2);
                    HFR_6km_dir_deg_going_to = mod(atan2(HFR_6km_E,HFR_6km_N)*180/pi+360,360);
                    if ~isnan(HFR_6km_speed_m_s)
                        HFR_6km_flag = true;
                    end
                end

                overall_flag = wind_flag & single_point_flag & ADCP_flag;% & HFR_2km_flag;

                if overall_flag

                    obs_wind = struct();
                    obs_current_profile = struct();
                    obs_wave_spectrum = struct();

                    % Wind
                    obs_wind.air_density_kg_m3 = air_density_kg_m3_value;
                    obs_wind.air_side_friction_velocity_m_s = ustar_m_s;
                    obs_wind.stress_angle_rad = stress_angle_rad;
                    obs_wind.wind_speed_m_s = windspeed_m_s;
                    obs_wind.wind_direction_rad = wind_dir_deg_coming_from*pi/180;
                    obs_wind.wind_meas_height_m = NDBC_struc.wind_height_m;

                    % Current_profile
                    z_m_obs = [ADCP_z(:); -1*single_point_depth_m];
                    current_speed_m_s_obs = [ADCP_speed_m_s(:); single_point_speed_m_s];
                    current_dir_rad_obs = [ADCP_dir_deg_going_to(:); single_point_dir_deg_going_to]*pi/180;
                    obs_current_profile.z_m = z_m_obs;
                    obs_current_profile.current_speed_m_s = current_speed_m_s_obs;
                    obs_current_profile.current_dir_rad = current_dir_rad_obs;

                    % Wave spectrum
                    f_Hz = f_Hz_vec(:);
                    theta_rad = theta_rad_vec(:)';
                    num_freq = length(f_Hz);
                    num_dir = length(theta_rad);
                    [s1,s2] = size(F_f_theta);
                    if s2 == num_freq
                        F_f_theta = F_f_theta';
                    end
                    F_f_theta_full = F_f_theta;
                    inds_retain = F_f_theta_full(:,1) > 0 & f_Hz < f_cut;
                    f_Hz = f_Hz(inds_retain);
                    F_f_theta_full = F_f_theta_full(inds_retain,:);
                    obs_wave_spectrum.f_Hz = f_Hz;
                    obs_wave_spectrum.theta_rad = theta_rad;
                    obs_wave_spectrum.F_f_theta = F_f_theta_full;

                    month_struc(spec_num).obs_wind = obs_wind;
                    month_struc(spec_num).obs_current_profile = obs_current_profile;
                    month_struc(spec_num).obs_wave_spectrum = obs_wave_spectrum;

                    if ~isnan(ustar_m_s)

                        % Perform computation with full current profile
                        fullcurrent_struc = compute_form_drag_effective_current_observation(obs_wind,obs_current_profile,obs_wave_spectrum,water_depth_m,k_max,true);

                        % Perform computation with only single point
                        % current meter data, 'slab' assumption
                        obs_current_profile.current_speed_m_s = obs_current_profile.current_speed_m_s*0 + single_point_speed_m_s;
                        obs_current_profile.current_dir_rad = obs_current_profile.current_dir_rad*0 + single_point_dir_deg_going_to*pi/180;
                        slabcurrent_struc = compute_form_drag_effective_current_observation(obs_wind,obs_current_profile,obs_wave_spectrum,water_depth_m,k_max,false);

                        % Perform computation with no current
                        obs_current_profile.current_speed_m_s = obs_current_profile.current_speed_m_s*0;
                        obs_current_profile.current_dir_rad = obs_current_profile.current_dir_rad*0;
                        nocurrent_struc = compute_form_drag_effective_current_observation(obs_wind,obs_current_profile,obs_wave_spectrum,water_depth_m,k_max,false);

                    else

                        downwind_stress_N_m2 = [];
                        crosswind_stress_N_m2 = [];

                    end

                    month_struc(spec_num).downwind_stress_N_m2 = downwind_stress_N_m2;
                    month_struc(spec_num).crosswind_stress_N_m2 = crosswind_stress_N_m2;

                    month_struc(spec_num).fullcurrent_struc = fullcurrent_struc;
                    month_struc(spec_num).slabcurrent_struc = slabcurrent_struc;
                    month_struc(spec_num).nocurrent_struc = nocurrent_struc;

                else

                    month_struc(spec_num).obs_wind = struct();
                    month_struc(spec_num).obs_current_profile = struct();
                    month_struc(spec_num).obs_wave_spectrum = struct();

                    month_struc(spec_num).downwind_stress_N_m2 = NaN;
                    month_struc(spec_num).crosswind_stress_N_m2 = NaN;

                    month_struc(spec_num).fullcurrent_struc = struct();
                    month_struc(spec_num).slabcurrent_struc = struct();
                    month_struc(spec_num).nocurrent_struc = struct();


                end

            else

                month_struc(spec_num).obs_wind = struct();
                month_struc(spec_num).obs_current_profile = struct();
                month_struc(spec_num).obs_wave_spectrum = struct();

                month_struc(spec_num).downwind_stress_N_m2 = NaN;
                month_struc(spec_num).crosswind_stress_N_m2 = NaN;

                month_struc(spec_num).fullcurrent_struc = struct();
                month_struc(spec_num).slabcurrent_struc = struct();
                month_struc(spec_num).nocurrent_struc = struct();

            end

        end

        s = struct();

        f_size = 70;
        k_size = 64;
        df = 0.005;

        f_NaN = NaN*ones(f_size,1);
        k_NaN = NaN*ones(k_size,1);

        fnames = {'DTime','air_density_kg_m3','wind_speed_m_s',...
            'wind_direction_rad','fp_Hz','Dp_rad','Dm_rad','f_Hz','F_f_m2_Hz','downwind_stress_N_m2',...
            'crosswind_stress_N_m2','k_rad_m_fullcurrent','F_k_m3_fullcurrent',...
            'tau_w_fullcurrent','tau_v_fullcurrent','tau_s_fullcurrent','U_w_fullcurrent',...
            'U_v_fullcurrent','U_drift_fullcurrent','k_rad_m_slabcurrent',...
            'F_k_m3_slabcurrent','tau_w_slabcurrent','tau_v_slabcurrent','tau_s_slabcurrent',...
            'U_w_slabcurrent','U_v_slabcurrent','U_drift_slabcurrent',...
            'k_rad_m_nocurrent','F_k_m3_nocurrent','tau_w_nocurrent',...
            'tau_v_nocurrent','tau_s_nocurrent','U_w_nocurrent','U_v_nocurrent','U_drift_nocurrent',...
            'current_z_m','current_U_m_s','current_D_deg'};

        for i = 1:length(fnames)
            s(1).(fnames{i}) = [];
        end

        for spec_num = 1:num_spectra

            if ~isempty(month_struc(spec_num).downwind_stress_N_m2) & ~isnan(month_struc(spec_num).downwind_stress_N_m2)

                s(spec_num).DTime = DTime_vec(spec_num);

                s(spec_num).air_density_kg_m3 = month_struc(spec_num).obs_wind.air_density_kg_m3;
                s(spec_num).wind_speed_m_s = month_struc(spec_num).obs_wind.wind_speed_m_s;
                s(spec_num).wind_direction_rad = month_struc(spec_num).obs_wind.wind_direction_rad;

                f_Hz = month_struc(spec_num).obs_wave_spectrum.f_Hz;
                theta_rad = month_struc(spec_num).obs_wave_spectrum.theta_rad;
                F_f_theta = month_struc(spec_num).obs_wave_spectrum.F_f_theta;
                F_f_theta_smooth = smoothdata2(F_f_theta,'movmean',{3,1},'omitnan');
                [r,c] = find(F_f_theta_smooth==max(F_f_theta_smooth,[],'all','omitnan'));
                fp = f_Hz(r);
                dp = theta_rad(c);
                dm = trapz(f_Hz,trapz(theta_rad,theta_rad.*F_f_theta,2))/trapz(f_Hz,trapz(theta_rad,F_f_theta,2));

                s(spec_num).fp_Hz = fp;
                s(spec_num).Dp_rad = dp;
                s(spec_num).Dm_rad = dm;

                f_bit = month_struc(spec_num).obs_wave_spectrum.f_Hz;
                F_f_bit = trapz(month_struc(spec_num).obs_wave_spectrum.theta_rad,month_struc(spec_num).obs_wave_spectrum.F_f_theta,2);

                L = length(f_bit);
                f_prepended = (df:df:df*(f_size-L))';

                s(spec_num).f_Hz = [f_prepended; f_bit];
                s(spec_num).F_f_m2_Hz = [0*f_prepended; F_f_bit];

                s(spec_num).downwind_stress_N_m2 = month_struc(spec_num).downwind_stress_N_m2;
                s(spec_num).crosswind_stress_N_m2 = month_struc(spec_num).crosswind_stress_N_m2;

                % Pull out full current profile with wind drift
                s(spec_num).current_z_m = month_struc(spec_num).fullcurrent_struc.current_z_m;
                s(spec_num).current_U_m_s = month_struc(spec_num).fullcurrent_struc.current_U_m_s;
                s(spec_num).current_D_deg = month_struc(spec_num).fullcurrent_struc.current_D_deg;

                % Full current profile, with wind drift
                s(spec_num).k_rad_m_fullcurrent = month_struc(spec_num).fullcurrent_struc.k_rad_m;
                s(spec_num).F_k_m3_fullcurrent = month_struc(spec_num).fullcurrent_struc.wave_F_k;

                s(spec_num).tau_w_fullcurrent = month_struc(spec_num).fullcurrent_struc.tau_w;
                s(spec_num).tau_v_fullcurrent = month_struc(spec_num).fullcurrent_struc.tau_v;
                s(spec_num).tau_s_fullcurrent = month_struc(spec_num).fullcurrent_struc.tau_s;

                s(spec_num).U_w_fullcurrent = month_struc(spec_num).fullcurrent_struc.U_w;
                s(spec_num).U_v_fullcurrent = month_struc(spec_num).fullcurrent_struc.U_v;

                s(spec_num).U_drift_fullcurrent = month_struc(spec_num).fullcurrent_struc.U_drift;

                % Slab current

                s(spec_num).k_rad_m_slabcurrent = month_struc(spec_num).slabcurrent_struc.k_rad_m;
                s(spec_num).F_k_m3_slabcurrent = month_struc(spec_num).slabcurrent_struc.wave_F_k;

                s(spec_num).tau_w_slabcurrent = month_struc(spec_num).slabcurrent_struc.tau_w;
                s(spec_num).tau_v_slabcurrent = month_struc(spec_num).slabcurrent_struc.tau_v;
                s(spec_num).tau_s_slabcurrent = month_struc(spec_num).slabcurrent_struc.tau_s;

                s(spec_num).U_w_slabcurrent = month_struc(spec_num).slabcurrent_struc.U_w;
                s(spec_num).U_v_slabcurrent = month_struc(spec_num).slabcurrent_struc.U_v;

                s(spec_num).U_drift_slabcurrent = month_struc(spec_num).slabcurrent_struc.U_drift;

                % No current
                s(spec_num).k_rad_m_nocurrent = month_struc(spec_num).nocurrent_struc.k_rad_m;
                s(spec_num).F_k_m3_nocurrent = month_struc(spec_num).nocurrent_struc.wave_F_k;

                s(spec_num).tau_w_nocurrent = month_struc(spec_num).nocurrent_struc.tau_w;
                s(spec_num).tau_v_nocurrent = month_struc(spec_num).nocurrent_struc.tau_v;
                s(spec_num).tau_s_nocurrent = month_struc(spec_num).nocurrent_struc.tau_s;

                s(spec_num).U_w_nocurrent = month_struc(spec_num).nocurrent_struc.U_w;
                s(spec_num).U_v_nocurrent = month_struc(spec_num).nocurrent_struc.U_v;

                s(spec_num).U_drift_nocurrent = month_struc(spec_num).nocurrent_struc.U_drift;

            end

        end

        month_struc = s;

        fnames = fieldnames(month_struc);

        clear s

        for i = 1:length(month_struc)

            if isempty(month_struc(i).DTime)

                for j = 1:length(fnames)

                    month_struc(i).(fnames{j}) = NaN;

                end

                month_struc(i).DTime = datetime(datevec(datenum(NaN)),'TimeZone','UTC');
                month_struc(i).f_Hz = NaN*ones(f_size,1);
                month_struc(i).F_f_m2_Hz = NaN*ones(f_size,1);

                month_struc(i).k_rad_m_fullcurrent = NaN*ones(k_size,1);
                month_struc(i).F_k_m3_fullcurrent = NaN*ones(k_size,1);

                month_struc(i).k_rad_m_slabcurrent = NaN*ones(k_size,1);
                month_struc(i).F_k_m3_slabcurrent = NaN*ones(k_size,1);

                month_struc(i).k_rad_m_nocurrent = NaN*ones(k_size,1);
                month_struc(i).F_k_m3_nocurrent = NaN*ones(k_size,1);

            end

        end

        DV_median = median(DV,'omitnan');
        YYYY = num2str(DV_median(1));
        MM = leading_zeroes(DV_median(2),2);

        out_name = [out_path YYYY '_' MM '_stress_decomposition_effective_current.mat'];
        save(out_name,'month_struc','-v7.3')

        clear month_struc

    end

    disp(['DONE WITH ' num2str(yearvec(year_ind))])

end
