%
function wind_stress_overestimation(fignum,fsize,numbins,quantiles)

load('_data/OOI_CE_data_compilation.mat')

EC_CD_w = EC_tau_w./air_density_kg_m3./(EC_U10_m_s-EC_U_w).^2;
EC_CD_v = EC_tau_v./air_density_kg_m3./(EC_U10_m_s-EC_U_v).^2;

label_x = 0.06;
label_y = 0.95;

labelcell = {'(a)','(b)','(c)','(d)'};

cmap = viridis(7);
violet = cmap(1,:);
teal = cmap(4,:);

msize = 10;

windlims = [0 12];

[EC_U10_m_s_fullcurrent_quantiles,EC_CD_w_fullcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(1,:),EC_CD_w(1,:),numbins,quantiles);
[EC_U10_m_s_slabcurrent_quantiles,EC_CD_w_slabcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(2,:),EC_CD_w(2,:),numbins,quantiles);
[EC_U10_m_s_nocurrent_quantiles,EC_CD_w_nocurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(3,:),EC_CD_w(3,:),numbins,quantiles);

EC_CD_w_fullcurrent_quantiles = 1000*EC_CD_w_fullcurrent_quantiles;
EC_CD_w_slabcurrent_quantiles = 1000*EC_CD_w_slabcurrent_quantiles;
EC_CD_w_nocurrent_quantiles = 1000*EC_CD_w_nocurrent_quantiles;

[~,EC_CD_v_fullcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(1,:),EC_CD_v(1,:),numbins,quantiles);
[~,EC_CD_v_slabcurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(2,:),EC_CD_v(2,:),numbins,quantiles);
[~,EC_CD_v_nocurrent_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s-EC_U_w(3,:),EC_CD_v(3,:),numbins,quantiles);

EC_CD_v_fullcurrent_quantiles = 1000*EC_CD_v_fullcurrent_quantiles;
EC_CD_v_slabcurrent_quantiles = 1000*EC_CD_v_slabcurrent_quantiles;
EC_CD_v_nocurrent_quantiles = 1000*EC_CD_v_nocurrent_quantiles;

U10_trial = EC_U10_m_s_nocurrent_quantiles(:,3);

rho_a_trial = 1.25;

tau_v_nocurrent = rho_a_trial*EC_CD_v_nocurrent_quantiles(:,3).*EC_U10_m_s_nocurrent_quantiles(:,3).^2/1000;
tau_v_slabcurrent = rho_a_trial*EC_CD_v_slabcurrent_quantiles(:,3).*EC_U10_m_s_slabcurrent_quantiles(:,3).^2/1000;
tau_v_fullcurrent = rho_a_trial*EC_CD_v_fullcurrent_quantiles(:,3).*EC_U10_m_s_fullcurrent_quantiles(:,3).^2/1000;

tau_w_nocurrent = rho_a_trial*EC_CD_w_nocurrent_quantiles(:,3).*EC_U10_m_s_nocurrent_quantiles(:,3).^2/1000;
tau_w_slabcurrent = rho_a_trial*EC_CD_w_slabcurrent_quantiles(:,3).*EC_U10_m_s_slabcurrent_quantiles(:,3).^2/1000;
tau_w_fullcurrent = rho_a_trial*EC_CD_w_fullcurrent_quantiles(:,3).*EC_U10_m_s_fullcurrent_quantiles(:,3).^2/1000;

tau_nocurrent = tau_v_nocurrent + tau_w_nocurrent;
tau_slabcurrent = tau_v_slabcurrent + tau_w_slabcurrent;
tau_fullcurrent = tau_v_fullcurrent + tau_w_fullcurrent;

figure(fignum);clf
tlayout = tiledlayout(1,3);
ax_struc = struct();

nexttile(1)
hold on
plot(windlims,0*windlims,'k--','linewidth',3)
plot(U10_trial,100*(tau_v_nocurrent-tau_v_slabcurrent)./tau_v_slabcurrent,'k-','linewidth',2);
plot(U10_trial,100*(tau_v_nocurrent-tau_v_slabcurrent)./tau_v_slabcurrent,'-','Color',violet,'linewidth',1);
h_slab = plot(U10_trial,100*(tau_v_nocurrent-tau_v_slabcurrent)./tau_v_slabcurrent,'s','markerfacecolor',violet,'markeredgecolor','k','markersize',msize);
plot(U10_trial,100*(tau_v_nocurrent-tau_v_fullcurrent)./tau_v_fullcurrent,'k-','linewidth',2);
plot(U10_trial,100*(tau_v_nocurrent-tau_v_fullcurrent)./tau_v_fullcurrent,'-','Color',teal,'linewidth',1);
h_full = plot(U10_trial,100*(tau_v_nocurrent-tau_v_fullcurrent)./tau_v_fullcurrent,'s','markerfacecolor',teal,'markeredgecolor','k','markersize',msize);
hold off
xlim(windlims)
ylim([-0.5 1]*50)
xlabel('U_{10} [m s^{-1}]')
ylabel('relative to no current [%]')
title('\tau_{\nu}')
ax_struc(1).ax = gca;

nexttile(2)
hold on
plot(windlims,0*windlims,'k--','linewidth',3)
plot(U10_trial,100*(tau_w_nocurrent-tau_w_slabcurrent)./tau_w_slabcurrent,'k-','linewidth',2);
plot(U10_trial,100*(tau_w_nocurrent-tau_w_slabcurrent)./tau_w_slabcurrent,'-','Color',violet,'linewidth',1);
h_slab = plot(U10_trial,100*(tau_w_nocurrent-tau_w_slabcurrent)./tau_w_slabcurrent,'s','markerfacecolor',violet,'markeredgecolor','k','markersize',msize);
plot(U10_trial,100*(tau_w_nocurrent-tau_w_fullcurrent)./tau_w_fullcurrent,'k-','linewidth',2);
plot(U10_trial,100*(tau_w_nocurrent-tau_w_fullcurrent)./tau_w_fullcurrent,'-','Color',teal,'linewidth',1);
h_full = plot(U10_trial,100*(tau_w_nocurrent-tau_w_fullcurrent)./tau_w_fullcurrent,'s','markerfacecolor',teal,'markeredgecolor','k','markersize',msize);
hold off
xlim(windlims)
ylim([-0.5 1]*50)
xlabel('U_{10} [m s^{-1}]')
ylabel('relative to no current [%]')
title('\tau_{w}')
ax_struc(2).ax = gca;

nexttile(3)
hold on
plot(windlims,0*windlims,'k--','linewidth',3)
plot(U10_trial,100*(tau_nocurrent-tau_slabcurrent)./tau_slabcurrent,'k-','linewidth',2);
plot(U10_trial,100*(tau_nocurrent-tau_slabcurrent)./tau_slabcurrent,'-','Color',violet,'linewidth',1);
h_slab = plot(U10_trial,100*(tau_nocurrent-tau_slabcurrent)./tau_slabcurrent,'s','markerfacecolor',violet,'markeredgecolor','k','markersize',msize);
plot(U10_trial,100*(tau_nocurrent-tau_fullcurrent)./tau_fullcurrent,'k-','linewidth',2);
plot(U10_trial,100*(tau_nocurrent-tau_fullcurrent)./tau_fullcurrent,'-','Color',teal,'linewidth',1);
h_full = plot(U10_trial,100*(tau_nocurrent-tau_fullcurrent)./tau_fullcurrent,'s','markerfacecolor',teal,'markeredgecolor','k','markersize',msize);
hold off
xlim(windlims)
ylim([-0.5 1]*50)
xlabel('U_{10} [m s^{-1}]')
ylabel('relative to no current [%]')
title('\tau_{total}')
H = [h_slab h_full];
L = {'slab current','sheared current'};
legend(H,L,'Location','Northeast')
ax_struc(3).ax = gca;

for i = 1:3
    ax_struc(i).ax.YTick = -50:25:50;
    ax_struc(i).ax.Box = 'on';
    nexttile(i)
    text(label_x*1.2,label_y,labelcell{i},'HorizontalAlignment','center','Units','normalized','FontSize',fsize)
end

tile_cleaner(ax_struc,tlayout)