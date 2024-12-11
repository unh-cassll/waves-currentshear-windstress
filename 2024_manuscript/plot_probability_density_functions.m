%
function plot_probability_density_functions(fignum,fsize)

load('_data/OOI_CE_data_compilation.mat')

rainbowcolors = flipud(spectral(7));
violet = rainbowcolors(1,:);
teal = rainbowcolors(2,:);

skycolors = sky(3);
darkblue = skycolors(3,:);

coppercolors = copper(3);
darkcopper = coppercolors(2,:);

fA = 0.2;

plims = [0 5];
Ulims = [0 15];
USlims = [0 1];
Tlims = [0 12];
Hlims = [0 6];

nbins = 100;

plot_x = 0.06;
plot_y = 0.96;

[counts_U10,centers_U10] = hist(EC_U10_m_s,nbins);

frac_U10 = counts_U10/sum(counts_U10);

cumu_vec = cumtrapz(centers_U10,counts_U10/trapz(centers_U10,counts_U10));
ind_10 = find(cumu_vec>0.1,1,'first');
ind_25 = find(cumu_vec>0.25,1,'first');
ind_50 = find(cumu_vec>0.5,1,'first');
ind_75 = find(cumu_vec>0.75,1,'first');
ind_90 = find(cumu_vec>0.9,1,'first');

figure(fignum);clf
tlayout = tiledlayout(2,2);
ax_struc = struct();

nexttile(1)
hold on
f_id = fill([centers_U10(ind_10:ind_90) fliplr(centers_U10(ind_10:ind_90))],100*[frac_U10(ind_10:ind_90) 0*frac_U10(ind_10:ind_90)],violet);
f_iq = fill([centers_U10(ind_25:ind_75) fliplr(centers_U10(ind_25:ind_75))],100*[frac_U10(ind_25:ind_75) 0*frac_U10(ind_25:ind_75)],violet);
plot(centers_U10(ind_50)*[1 1],[0 100*frac_U10(ind_50)],'-','Color',violet,'LineWidth',3)
plot(centers_U10,frac_U10*100,'k-','LineWidth',3)
hold off
xlim(Ulims)
ylim(plims)
xlabel('U_{10} [m s^{-1}]')
text(plot_x,plot_y,'(a)','FontSize',fsize,'HorizontalAlignment','center','Units','normalized')

f_id.FaceAlpha = fA;
f_iq.FaceAlpha = fA;

ax_struc(1).ax = gca;

[counts_Us,centers_Us] = hist(U_sp_m_s,nbins);

frac_Us = counts_Us/sum(counts_Us);

cumu_vec = cumtrapz(centers_Us,counts_Us/trapz(centers_Us,counts_Us));
ind_10 = find(cumu_vec>0.1,1,'first');
ind_25 = find(cumu_vec>0.25,1,'first');
ind_50 = find(cumu_vec>0.5,1,'first');
ind_75 = find(cumu_vec>0.75,1,'first');
ind_90 = find(cumu_vec>0.9,1,'first');

nexttile(2)
hold on
f_id = fill([centers_Us(ind_10:ind_90) fliplr(centers_Us(ind_10:ind_90))],100*[frac_Us(ind_10:ind_90) 0*frac_Us(ind_10:ind_90)],darkblue);
f_iq = fill([centers_Us(ind_25:ind_75) fliplr(centers_Us(ind_25:ind_75))],100*[frac_Us(ind_25:ind_75) 0*frac_Us(ind_25:ind_75)],darkblue);
plot(centers_Us(ind_50)*[1 1],[0 100*frac_Us(ind_50)],'-','Color',darkblue,'LineWidth',3)
plot(centers_Us,frac_Us*100,'k-','LineWidth',3)
hold off
xlim(USlims)
ylim(plims)
xlabel('U_S [m s^{-1}]')
text(plot_x,plot_y,'(b)','FontSize',fsize,'HorizontalAlignment','center','Units','normalized')

f_id.FaceAlpha = fA;
f_iq.FaceAlpha = fA;

ax_struc(2).ax = gca;

[counts_Hs,centers_Hs] = hist(Hs,nbins);

frac_Hs = counts_Hs/sum(counts_Hs);

cumu_vec = cumtrapz(centers_Hs,counts_Hs/trapz(centers_Hs,counts_Hs));
ind_10 = find(cumu_vec>0.1,1,'first');
ind_25 = find(cumu_vec>0.25,1,'first');
ind_50 = find(cumu_vec>0.5,1,'first');
ind_75 = find(cumu_vec>0.75,1,'first');
ind_90 = find(cumu_vec>0.9,1,'first');

nexttile(3)
hold on
f_id = fill([centers_Hs(ind_10:ind_90) fliplr(centers_Hs(ind_10:ind_90))],100*[frac_Hs(ind_10:ind_90) 0*frac_Hs(ind_10:ind_90)],darkcopper);
f_iq = fill([centers_Hs(ind_25:ind_75) fliplr(centers_Hs(ind_25:ind_75))],100*[frac_Hs(ind_25:ind_75) 0*frac_Hs(ind_25:ind_75)],darkcopper);
plot(centers_Hs(ind_50)*[1 1],[0 100*frac_Hs(ind_50)],'-','Color',darkcopper,'LineWidth',3)
plot(centers_Hs,frac_Hs*100,'k-','LineWidth',3)
hold off
xlim(Hlims)
ylim(plims)
xlabel('H_s [m]')
text(plot_x,plot_y,'(c)','FontSize',fsize,'HorizontalAlignment','center','Units','normalized')

f_id.FaceAlpha = fA;
f_iq.FaceAlpha = fA;

ax_struc(3).ax = gca;

[counts_Tm,centers_Tm] = hist(Tm,nbins);

frac_Tm = counts_Tm/sum(counts_Tm);

cumu_vec = cumtrapz(centers_Tm,counts_Tm/trapz(centers_Tm,counts_Tm));
ind_10 = find(cumu_vec>0.1,1,'first');
ind_25 = find(cumu_vec>0.25,1,'first');
ind_50 = find(cumu_vec>0.5,1,'first');
ind_75 = find(cumu_vec>0.75,1,'first');
ind_90 = find(cumu_vec>0.9,1,'first');

nexttile(4)
hold on
f_id = fill([centers_Tm(ind_10:ind_90) fliplr(centers_Tm(ind_10:ind_90))],100*[frac_Tm(ind_10:ind_90) 0*frac_Tm(ind_10:ind_90)],teal);
f_iq = fill([centers_Tm(ind_25:ind_75) fliplr(centers_Tm(ind_25:ind_75))],100*[frac_Tm(ind_25:ind_75) 0*frac_Tm(ind_25:ind_75)],teal);
plot(centers_Tm(ind_50)*[1 1],[0 100*frac_Tm(ind_50)],'-','Color',teal,'LineWidth',3)
plot(centers_Tm,frac_Tm*100,'k-','LineWidth',3)
hold off
xlim(Tlims)
ylim(plims)
xlabel('T_m [s]')
text(plot_x,plot_y,'(d)','FontSize',fsize,'HorizontalAlignment','center','Units','normalized')

f_id.FaceAlpha = fA;
f_iq.FaceAlpha = fA;

ax_struc(4).ax = gca;

for i = 1:4
    ax_struc(i).ax.Box = 'on';
    if ~mod(i,2)
        ax_struc(i).ax.YTickLabel = '';
    else
        nexttile(i)
        ylabel('% of total observations')
    end
end

tlayout.TileSpacing = 'compact';