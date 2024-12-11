%
function stress_fraction(fignum,fsize,numbins,quantiles)

load('_data/OOI_CE_data_compilation.mat')

EC_tau_v_fraction = EC_tau_v./EC_stress_N_m2;
EC_tau_w_fraction = EC_tau_w./EC_stress_N_m2;
EC_tau_s_fraction = EC_tau_s./EC_stress_N_m2;

cmap = viridis(7);
violet = cmap(1,:);
teal = cmap(4,:);
crimson = [0.7 0 0];

lw = 3;

windlims = [0 12];

[EC_U10_m_s_fullcurrent_quantiles,EC_tau_w_fraction_fullcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(1,:),EC_tau_w_fraction(1,:),numbins,quantiles);
[EC_U10_m_s_slabcurrent_quantiles,EC_tau_w_fraction_slabcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(2,:),EC_tau_w_fraction(2,:),numbins,quantiles);
[EC_U10_m_s_nocurrent_quantiles,EC_tau_w_fraction_nocurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(3,:),EC_tau_w_fraction(3,:),numbins,quantiles);

[~,EC_tau_v_fraction_fullcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(1,:),EC_tau_v_fraction(1,:),numbins,quantiles);
[~,EC_tau_v_fraction_slabcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(2,:),EC_tau_v_fraction(2,:),numbins,quantiles);
[~,EC_tau_v_fraction_nocurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(3,:),EC_tau_v_fraction(3,:),numbins,quantiles);

[~,EC_tau_s_fraction_fullcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(1,:),EC_tau_s_fraction(1,:),numbins,quantiles);
[~,EC_tau_s_fraction_slabcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(2,:),EC_tau_s_fraction(2,:),numbins,quantiles);
[~,EC_tau_s_fraction_nocurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(3,:),EC_tau_s_fraction(3,:),numbins,quantiles);

figure(fignum);clf
hold on
h_v = plot([0 1],[-1 -1],'k-','linewidth',lw);
h_w = plot([0 1],[-1 -1],'k--','linewidth',lw);
h_s = plot([0 1],[-1 -1],'k:','linewidth',lw);
plot(EC_U10_m_s_nocurrent_quantiles(:,3),EC_tau_v_fraction_nocurrent_quantiles(:,3),'-','Color',crimson,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_slabcurrent_quantiles(:,3),EC_tau_v_fraction_slabcurrent_quantiles(:,3),'-','Color',violet,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_fullcurrent_quantiles(:,3),EC_tau_v_fraction_fullcurrent_quantiles(:,3),'-','Color',teal,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_nocurrent_quantiles(:,3),EC_tau_w_fraction_nocurrent_quantiles(:,3),'--','Color',crimson,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_slabcurrent_quantiles(:,3),EC_tau_w_fraction_slabcurrent_quantiles(:,3),'--','Color',violet,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_fullcurrent_quantiles(:,3),EC_tau_w_fraction_fullcurrent_quantiles(:,3),'--','Color',teal,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_nocurrent_quantiles(:,3),EC_tau_s_fraction_nocurrent_quantiles(:,3),':','Color',crimson,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_slabcurrent_quantiles(:,3),EC_tau_s_fraction_slabcurrent_quantiles(:,3),':','Color',violet,'markersize',15,'linewidth',lw);
plot(EC_U10_m_s_fullcurrent_quantiles(:,3),EC_tau_s_fraction_fullcurrent_quantiles(:,3),':','Color',teal,'markersize',15,'linewidth',lw);
hold off
box on
H = [h_v h_w h_s];
L = {'\tau_\nu','\tau_w','\tau_b'};
legend(H,L,'Location','northeast')
xlim(windlims)
ylim([0 1])
xlabel('U-U_{i} [m s^{-1}]')
ylabel('\tau_i/\tau')
text(0.6,0.57,'no current','Color',crimson,'FontSize',fsize,'HorizontalAlignment','left')
text(0.6,0.5,'slab current','Color',violet,'FontSize',fsize,'HorizontalAlignment','left')
text(0.6,0.43,'sheared current','Color',teal,'FontSize',fsize,'HorizontalAlignment','left')