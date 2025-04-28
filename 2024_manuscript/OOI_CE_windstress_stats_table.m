% Compute t-tests

close all;clear;clc

load('OOI_CE_windstress_data.mat')

CD_w = tau_w_all./rho_a_all./(U10_m_s_all-U_w_all).^2;
CD_v = tau_v_all./rho_a_all./(U10_m_s_all-U_v_all).^2;

CD_w_nocurrent = CD_w(3,:);
CD_v_nocurrent = CD_v(3,:);

U_lims = 0.5:2:13;

tau_v_nocurrent = rho_a_all.*CD_v_nocurrent.*(U10_m_s_all-U_v_all(3,:)).^2;
tau_w_nocurrent = rho_a_all.*CD_w_nocurrent.*(U10_m_s_all-U_w_all(3,:)).^2;

tau_v_slabcurrent = rho_a_all.*CD_v_nocurrent.*(U10_m_s_all-U_v_all(2,:)).^2;
tau_w_slabcurrent = rho_a_all.*CD_w_nocurrent.*(U10_m_s_all-U_w_all(2,:)).^2;

tau_v_fullcurrent = rho_a_all.*CD_v_nocurrent.*(U10_m_s_all-U_v_all(1,:)).^2;
tau_w_fullcurrent = rho_a_all.*CD_w_nocurrent.*(U10_m_s_all-U_w_all(1,:)).^2;

tau_nocurrent_block = real([tau_v_nocurrent; tau_w_nocurrent; tau_v_nocurrent+tau_w_nocurrent]);
tau_slabcurrent_block = real([tau_v_slabcurrent; tau_w_slabcurrent; tau_v_slabcurrent+tau_w_slabcurrent]);
tau_fullcurrent_block = real([tau_v_fullcurrent; tau_w_fullcurrent; tau_v_fullcurrent+tau_w_fullcurrent]);

s = struct();

for i = 1:length(U_lims)-1

    U_low = U_lims(i);
    U_high = U_lims(i+1);

    inds = find(U10_m_s_all > U_low & U10_m_s_all <= U_high & ~isnan(mean(tau_nocurrent_block+tau_slabcurrent_block+tau_fullcurrent_block)));

    U_snip = U10_m_s_all(inds);
    tau_nocurrent_snip = tau_nocurrent_block(:,inds);
    tau_slabcurrent_snip = tau_slabcurrent_block(:,inds);
    tau_fullcurrent_snip = tau_fullcurrent_block(:,inds);

    s(i).U_mean = mean(U_snip);
    s(i).N = length(U_snip);
	pd_full_no = pdiff(tau_nocurrent_snip(3,:),tau_fullcurrent_snip(3,:));
	pd_full_slab = pdiff(tau_slabcurrent_snip(3,:),tau_fullcurrent_snip(3,:));
	pd_slab_no = pdiff(tau_nocurrent_snip(3,:),tau_slabcurrent_snip(3,:));
	s(i).median_rel_diff = [pd_full_no pd_full_slab pd_slab_no]';
	s(i).tauhat = [...
	mean(tau_nocurrent_snip(3,:))/mean(tau_fullcurrent_snip(3,:)),...
	mean(tau_slabcurrent_snip(3,:))/mean(tau_fullcurrent_snip(3,:)),...
	mean(tau_fullcurrent_snip(3,:))/mean(tau_fullcurrent_snip(3,:))]';
	
    % Pair #1: sheared vs. no current
    [~,p_full_no] = ttest2(tau_fullcurrent_snip(3,:),tau_nocurrent_snip(3,:),'Vartype','unequal');

    % Pair #2: sheared vs. slab
    [~,p_full_slab] = ttest2(tau_fullcurrent_snip(3,:),tau_slabcurrent_snip(3,:),'Vartype','unequal');

    % Pair #3: slab vs. no
    [~,p_slab_no] = ttest2(tau_slabcurrent_snip(3,:),tau_nocurrent_snip(3,:),'Vartype','unequal');

    s(i).p = [p_full_no p_full_slab p_slab_no]';

end

bin_number = (1:length(U_lims)-1)';
U10_m_s_binned = [s.U_mean]';
p = [s.p]';
N = [s.N]';
reldiff = [s.median_rel_diff]';
tauhat = [s.tauhat]';

variable_names = {'bin number','U_10 [m/s]','p_full_no','p_full_slab','p_slab_no','N'};
table_of_ttest_results = table(bin_number,U10_m_s_binned,p(:,1),p(:,2),p(:,3),N,'VariableNames',variable_names)

variable_names = {'bin number','U_10 [m/s]','reldiff_full_no','reldiff_full_slab','reldiff_slab_no','N'};
table_of_reldiff_results = table(bin_number,U10_m_s_binned,reldiff(:,1),reldiff(:,2),reldiff(:,3),N,'VariableNames',variable_names)

variable_names = {'bin number','U_10 [m/s]','tauhat_full_no','tauhat_full_slab','tauhat_full_full','N'};
table_of_tauhat_results = table(bin_number,U10_m_s_binned,tauhat(:,1),tauhat(:,2),tauhat(:,3),N,'VariableNames',variable_names)
% 
function pd = pdiff(X,Y)
	pd = (mean(X)-mean(Y))/mean(Y);
end