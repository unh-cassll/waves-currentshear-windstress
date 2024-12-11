%
function drag_coefficient_variation(fignum,fsize,numbins,quantiles)

load('_data/OOI_CE_data_compilation.mat')

EC_CD_w = EC_tau_w./air_density_kg_m3./(EC_U10_m_s-EC_U_w).^2;
EC_CD_v = EC_tau_v./air_density_kg_m3./(EC_U10_m_s-EC_U_v).^2;

label_x = 0.06;
label_y = 0.95;

labelcell = {'(a)','(b)','(c)','(d)'};

lw_interquartile = 4;

cmap = viridis(7);
violet = cmap(1,:);
teal = cmap(4,:);
crimson = [0.7 0 0];

windlims = [0 12];
CD_var_lims = [1 100];
Cdlims = [1e-2 1e1]*2;

CD_var_ticks = [1 10 100];
CD_var_ticklabels = {'1','10','100'};

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

figure(fignum);clf
tlayout = tiledlayout(2,2);
nexttile(1)
hold on
f_nocurrent = fill([EC_U10_m_s_nocurrent_quantiles(:,3); flipud(EC_U10_m_s_nocurrent_quantiles(:,3))],[EC_CD_v_nocurrent_quantiles(:,2); flipud(EC_CD_v_nocurrent_quantiles(:,4))],crimson);
h_nocurrent = plot(EC_U10_m_s_nocurrent_quantiles(:,3),EC_CD_v_nocurrent_quantiles(:,3),'-','Color',crimson,'linewidth',lw_interquartile);
f_slabcurrent = fill([EC_U10_m_s_slabcurrent_quantiles(:,3); flipud(EC_U10_m_s_slabcurrent_quantiles(:,3))],[EC_CD_v_slabcurrent_quantiles(:,2); flipud(EC_CD_v_slabcurrent_quantiles(:,4))],violet);
h_slabcurrent = plot(EC_U10_m_s_slabcurrent_quantiles(:,3),EC_CD_v_slabcurrent_quantiles(:,3),'-','Color',violet,'linewidth',lw_interquartile);
f_fullcurrent = fill([EC_U10_m_s_fullcurrent_quantiles(:,3); flipud(EC_U10_m_s_fullcurrent_quantiles(:,3))],[EC_CD_v_fullcurrent_quantiles(:,2); flipud(EC_CD_v_fullcurrent_quantiles(:,4))],teal);
h_fullcurrent = plot(EC_U10_m_s_fullcurrent_quantiles(:,3),EC_CD_v_fullcurrent_quantiles(:,3),'-','Color',teal,'linewidth',lw_interquartile);
hold off
box on
xlim(windlims)
ylim(Cdlims)
xlabel('U-U_{i} [m s^{-1}]')
ylabel('1000\timesC_{D,i}')
f_nocurrent.FaceAlpha = 0.2;
f_nocurrent.LineStyle = 'none';
f_slabcurrent.FaceAlpha = 0.2;
f_slabcurrent.LineStyle = 'none';
f_fullcurrent.FaceAlpha = 0.2;
f_fullcurrent.LineStyle = 'none';
ax_struc(1).ax = gca;
ax_struc(1).ax.YScale = 'log';
title('C_{D,\nu}')

nexttile(2)
hold on
f_nocurrent = fill([EC_U10_m_s_nocurrent_quantiles(:,3); flipud(EC_U10_m_s_nocurrent_quantiles(:,3))],[EC_CD_w_nocurrent_quantiles(:,2); flipud(EC_CD_w_nocurrent_quantiles(:,4))],crimson);
h_nocurrent = plot(EC_U10_m_s_nocurrent_quantiles(:,3),EC_CD_w_nocurrent_quantiles(:,3),'-','Color',crimson,'linewidth',lw_interquartile);
f_slabcurrent = fill([EC_U10_m_s_slabcurrent_quantiles(:,3); flipud(EC_U10_m_s_slabcurrent_quantiles(:,3))],[EC_CD_w_slabcurrent_quantiles(:,2); flipud(EC_CD_w_slabcurrent_quantiles(:,4))],violet);
h_slabcurrent = plot(EC_U10_m_s_slabcurrent_quantiles(:,3),EC_CD_w_slabcurrent_quantiles(:,3),'-','Color',violet,'linewidth',lw_interquartile);
f_fullcurrent = fill([EC_U10_m_s_fullcurrent_quantiles(:,3); flipud(EC_U10_m_s_fullcurrent_quantiles(:,3))],[EC_CD_w_fullcurrent_quantiles(:,2); flipud(EC_CD_w_fullcurrent_quantiles(:,4))],teal);
h_fullcurrent = plot(EC_U10_m_s_fullcurrent_quantiles(:,3),EC_CD_w_fullcurrent_quantiles(:,3),'-','Color',teal,'linewidth',lw_interquartile);
hold off
box on
xlim(windlims)
ylim(Cdlims)
xlabel('U-U_{i} [m s^{-1}]')
ylabel('1000\timesC_{D,i}')
H = [h_nocurrent h_slabcurrent h_fullcurrent];
L = {'no current','slab current','sheared current'};
legend(H,L,'location','northeast')
f_nocurrent.FaceAlpha = 0.2;
f_nocurrent.LineStyle = 'none';
f_slabcurrent.FaceAlpha = 0.2;
f_slabcurrent.LineStyle = 'none';
f_fullcurrent.FaceAlpha = 0.2;
f_fullcurrent.LineStyle = 'none';
ax_struc(2).ax = gca;
ax_struc(2).ax.YScale = 'log';
title('C_{D,w}')

for i = 1:2
    ax_struc(i).ax.YTick = 10.^(-2:1:1);
    ax_struc(i).ax.YTickLabel = {'0.01','0.1','1','10'};
end

nexttile(3)
hold on
h_nocurrent = plot(EC_U10_m_s_nocurrent_quantiles(:,3),10.^(log10(EC_CD_v_nocurrent_quantiles(:,4))-log10(EC_CD_v_nocurrent_quantiles(:,2))),'-','Color',crimson,'linewidth',lw_interquartile);
h_slabcurrent = plot(EC_U10_m_s_slabcurrent_quantiles(:,3),10.^(log10(EC_CD_v_slabcurrent_quantiles(:,4))-log10(EC_CD_v_slabcurrent_quantiles(:,2))),'-','Color',violet,'linewidth',lw_interquartile);
h_fullcurrent = plot(EC_U10_m_s_fullcurrent_quantiles(:,3),10.^(log10(EC_CD_v_fullcurrent_quantiles(:,4))-log10(EC_CD_v_fullcurrent_quantiles(:,2))),'-','Color',teal,'linewidth',lw_interquartile);
hold off
box on
xlim(windlims)
ylim(CD_var_lims)
ax_struc(3).ax=gca;
ax_struc(3).ax.YScale = 'log';
ax_struc(3).ax.YTick = CD_var_ticks;
ax_struc(3).ax.YTickLabels = CD_var_ticklabels;
xlabel('U-U_{i} [m s^{-1}]')
ylabel('C_{D,i} interquartile range size')

nexttile(4)
hold on
h_nocurrent = plot(EC_U10_m_s_nocurrent_quantiles(:,3),10.^(log10(EC_CD_w_nocurrent_quantiles(:,4))-log10(EC_CD_w_nocurrent_quantiles(:,2))),'-','Color',crimson,'linewidth',lw_interquartile);
h_slabcurrent = plot(EC_U10_m_s_slabcurrent_quantiles(:,3),10.^(log10(EC_CD_w_slabcurrent_quantiles(:,4))-log10(EC_CD_w_slabcurrent_quantiles(:,2))),'-','Color',violet,'linewidth',lw_interquartile);
h_fullcurrent = plot(EC_U10_m_s_fullcurrent_quantiles(:,3),10.^(log10(EC_CD_w_fullcurrent_quantiles(:,4))-log10(EC_CD_w_fullcurrent_quantiles(:,2))),'-','Color',teal,'linewidth',lw_interquartile);
hold off
box on
xlim(windlims)
ylim(CD_var_lims)
ax_struc(4).ax=gca;
ax_struc(4).ax.YScale = 'log';
ax_struc(4).ax.YTick = CD_var_ticks;
ax_struc(4).ax.YTickLabels = CD_var_ticklabels;
xlabel('U-U_{i} [m s^{-1}]')
for i = 1:4
    nexttile(i)
    text(label_x,label_y,labelcell{i},'HorizontalAlignment','center','Units','normalized','FontSize',fsize)
end

tile_cleaner(ax_struc,tlayout)