%
function make_maps_PNW_CE_OOI(fignum,fsize)

load('_data/CE_OOI_bathy.mat')

lon_center = -124.304;
lat_center = 44.6393;

latlims = [43.75 48.5];
lonlims = [-128 -121];

patchcolor = [253 190 110]/255;

rusty = [0.75 0 0];

lonlims_inset = [-124.5 -123.9167];
latlims_inset = [44.4 44.8];

dz = 20;
depth_vals = -120:dz:0;
depth_vals_str = num2cell(depth_vals);
n_levels = length(depth_vals);

cmap = deep(n_levels);

mask = NaN*depth_m;
mask(depth_m < 0) = 1;

figure(fignum);clf
tlayout = tiledlayout(1,2);

nexttile(1)
m_proj('miller','long',lonlims,'lat',latlims);
m_gshhs_f('patch',patchcolor);
m_gshhs_f('speckle','color','k');
hold on
m_plot([lonlims_inset(1) lonlims_inset(2) lonlims_inset(2) lonlims_inset(1) lonlims_inset(1)],[latlims_inset(1) latlims_inset(1) latlims_inset(2) latlims_inset(2) latlims_inset(1)],'-','Color',rusty,'LineWidth',2)
m_plot(lon_center,lat_center,'o','markeredgecolor','k','markerfacecolor',rusty,'markersize',5)
hold off
m_text(-127.75,48.25,'(a)','FontSize',fsize)
m_grid('box','fancy','fontsize',fsize,'linewidth',2,'tickdir','out','xaxisloc','bottom','yaxisloc','left');

nexttile(2)
m_proj('miller','long',lonlims_inset,'lat',latlims_inset);
m_gshhs_f('patch',patchcolor);
m_gshhs_f('speckle','color','k');
hold on
m_contour(lon_vec,lat_vec,(mask.*depth_m)',depth_vals,'LineWidth',2);colormap(cmap(1:end-1,:))
m_plot(lon_center,lat_center,'o','markeredgecolor','k','markerfacecolor',rusty,'markersize',7)
hold off
m_text(-124.48,44.775,'(b)','FontSize',fsize)
m_grid('box','fancy','fontsize',fsize,'linewidth',2,'tickdir','out','xaxisloc','bottom','yaxisloc','right');
cbar = colorbar;
cbar.Ticks = depth_vals;
cbar.TickLabels = depth_vals_str(1:end-1);
clim([depth_vals(1) depth_vals(end-1)]+[-1 1]*dz/2)
set(get(cbar,'Title'),'String','depth [m]')
cbar.Location = 'northoutside';