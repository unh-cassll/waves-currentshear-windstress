%
function plot_roses(fignum,fsize)

load('_data/OOI_CE_data_compilation.mat')

Options = struct();
Options.AngleNorth = 0;
Options.AngleEast = 90;
Options.axesfontname = 'Liberation Serif';
Options.textfontname = 'Liberation Serif';
Options.frequencyFontSize = fsize*3/4;
Options.axesfontsize = fsize;
Options.titlefontsize = fsize;
Options.legendfontsize = fsize*3/4;
Options.TitleString = '';
Options.legendvariable = 'U_{10}';
Options.lablegend = [];
Options.legendposition = 'southeast';
Options.FreqLabelAngle = 45;
Options.NSpeeds = 6;

U_high = 15;
inds_keep = EC_U10_m_s <= U_high;

Options.cmap = flipud(spectral);

figure(fignum);clf
tlayout = tiledlayout(2,2);

nexttile()
ax_struc(1).ax = gca;
Options.axes = ax_struc(1).ax;
WindRose(mod(wind_direction_rad(inds_keep)*180/pi,360),EC_U10_m_s(inds_keep),Options);

Options.cmap = deep;

Options.legendvariable = 'U_S';

sp_high = 2;
inds_keep = U_sp_m_s <= sp_high;

nexttile()
ax_struc(2).ax = gca;
Options.axes = ax_struc(2).ax;
WindRose(D_sp_rad(inds_keep)*180/pi+180,U_sp_m_s(inds_keep),Options);

Options.cmap = copper;

Options.legendvariable = 'H_s';

Hs_high = 6;
inds_keep = Hs <= Hs_high;

nexttile()
ax_struc(3).ax = gca;
Options.axes = ax_struc(3).ax;
WindRose(Dm_rad(inds_keep)*180/pi+180,Hs(inds_keep),Options);

Options.cmap = viridis;

Options.legendvariable = 'T_m';

Tm_high = 10;
inds_keep = Tm <= Tm_high;

nexttile()
ax_struc(4).ax = gca;
Options.axes = ax_struc(4).ax;
WindRose(Dm_rad(inds_keep)*180/pi+180,Tm(inds_keep),Options);
