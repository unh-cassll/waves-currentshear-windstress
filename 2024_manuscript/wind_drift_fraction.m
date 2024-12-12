%
function wind_drift_fraction(fignum,fsize,numbins,quantiles)

load('_data/OOI_CE_data_compilation.mat')

label_x = 0.06;
label_y = 0.95;

labelcell = {'(a)','(b)','(c)','(d)'};

lw_interdecile = 2;

cmap = viridis(7);
bluish = cmap(3,:);
windlims = [0 12];

[EC_U10_m_s_nocurrent_quantiles,EC_U_drift_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s,EC_U_drift,numbins,quantiles);
[~,EC_ustar_quantiles] = compute_quantiles_fixed_binsize(EC_U10_m_s,EC_air_friction_velocity_m_s,numbins,quantiles);

figure(fignum);clf

tlayout = tiledlayout(2,1);
ax_struc = struct();
nexttile(1)
hold on
f_interquartile = fill([EC_U10_m_s_nocurrent_quantiles(:,3); flipud(EC_U10_m_s_nocurrent_quantiles(:,3))],100*[EC_U_drift_quantiles(:,2)./EC_U10_m_s_nocurrent_quantiles(:,3); flipud(EC_U_drift_quantiles(:,4)./EC_U10_m_s_nocurrent_quantiles(:,3))],bluish);
plot(EC_U10_m_s_nocurrent_quantiles(:,3),100*EC_U_drift_quantiles(:,3)./EC_U10_m_s_nocurrent_quantiles(:,3),'-','Color',bluish,'linewidth',lw_interdecile)
plot(windlims,0*windlims+3,'k--','linewidth',2)
hold off
box on
xlabel('$U_{10}\,[m\,s^{-1}]$','Interpreter','LaTeX')
ylabel('$100\times\mathcal{U}_{S}/U_{10}\,[\%]$','Interpreter','LaTeX')
text(0.86*windlims(end),2.3,'$3\%\times U_{10}$','HorizontalAlignment','center','FontSize',fsize,'Interpreter','LaTeX')
xlim(windlims);ylim([0 10])
f_interquartile.LineStyle = 'none';
f_interquartile.FaceAlpha = 0.2;
ax_struc(1).ax = gca;

nexttile(2)
hold on
f_interquartile = fill([EC_U10_m_s_nocurrent_quantiles(:,3); flipud(EC_U10_m_s_nocurrent_quantiles(:,3))],[EC_U_drift_quantiles(:,2)./EC_ustar_quantiles(:,3); flipud(EC_U_drift_quantiles(:,4)./EC_ustar_quantiles(:,3))],bluish);
plot(EC_U10_m_s_nocurrent_quantiles(:,3),EC_U_drift_quantiles(:,3)./EC_ustar_quantiles(:,3),'-','Color',bluish,'linewidth',lw_interdecile)
plot(windlims,0*windlims+0.53,'k--','linewidth',2)
hold off
box on
xlabel('$U_{10}\,[m\,s^{-1}]$','Interpreter','LaTeX')
ylabel('$\mathcal{U}_{S}/u_{*}$','Interpreter','LaTeX')
text(0.9*windlims(end),0.45,'$0.53u_*$','HorizontalAlignment','center','FontSize',fsize,'Interpreter','LaTeX')
xlim(windlims);ylim([0 1.4])
f_interquartile.LineStyle = 'none';
f_interquartile.FaceAlpha = 0.2;
ax_struc(2).ax = gca;

for i = 1:length(ax_struc)
    nexttile(i)
    text(label_x,label_y,labelcell{i},'HorizontalAlignment','center','Units','normalized','FontSize',fsize)
end

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';