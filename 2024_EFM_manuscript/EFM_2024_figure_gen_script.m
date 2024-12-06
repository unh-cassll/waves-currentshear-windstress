% Figure generation script for
% "Accounting for Waves and Current Shear in Ocean Wind Stress Parameterization"
% 
% D. G. Ortiz-Suslow, N. J. M. Laxague, M. Curcic, & J.-V. I. Bjorkqvist, 2024
%

% Figure 3 uses the M_MAP toolbox, which you can acquire here:
% https://www.eoas.ubc.ca/~rich/map.html
% Once you've downloaded the ZIP, add 'm_map' to your path

addpath _codes
addpath _data
addpath _outputs

% addpath ../../../../geographic_utilities/
addpath ../../../../geographic_utilities/m_map/

close all;clear;clc

fsize = 22;

figure_style(fsize)

numbins = 12;
quantiles = [10 25 50 75 90];

n = 1;
print_options = {'none','svg','png'};
print_option = print_options{n};

dpi_val = 100;
dpi_string = ['-r' num2str(dpi_val)];

boxpos = [1250 550 600 600];
% two_by_one_pos = [25 25 600 1200];
two_by_one_pos = [1250 550 600 1200];

out_path = '_outputs/';

%% Figure 1 - conceptual diagram (drawn in InkScape)

%% Figure 2 - wind-induced shear profiles

%% Figure 3 - map of PNW and CE OOI site

fignum = 3;
figure(fignum)
set(gcf,'Position',boxpos.*[1 1 2 1].*[1 1 1.5 1.5])

make_maps_PNW_CE_OOI(fignum,fsize)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'PNW_CE_maps.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'PNW_CE_maps.png'],'-dpng',dpi_string);
end

%% Figure 4 - "roses" (direction/magnitude histograms)

fignum = 4;
figure(fignum)
set(gcf,'Position',boxpos.*[1 1 2 2])

plot_roses(fignum,fsize)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'bouquet_of_roses.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'bouquet_of_roses.png'],'-dpng',dpi_string);
end

%% Figure 5 - probability density functions

fignum = 5;
figure(fignum)
set(gcf,'Position',boxpos.*[1 1 2 2])

plot_probability_density_functions(fignum,fsize)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'histogram_quartet.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'histogram_quartet.png'],'-dpng',dpi_string);
end

%% Figure 6 - wind drift fraction

fignum = 6;
figure(fignum)
set(gcf,'Position',two_by_one_pos)

wind_drift_fraction(fignum,fsize,numbins,quantiles)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'Udrift_fraction.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'Udrift_fraction.png'],'-dpng',dpi_string);
end

%% Figure 7 - stress decomposition

fignum = 7;
figure(fignum)
set(gcf,'Position',boxpos)

stress_fraction(fignum,fsize,numbins,quantiles)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'wave_stress_fraction_MV2009.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'wave_stress_fraction_MV2009.png'],'-dpng',dpi_string);
end


%% Figure 8 - variation of C_D: no current, slab current, sheared current

fignum = 8;
figure(fignum)
set(gcf,'Position',boxpos.*[1 1 2 2])

drag_coefficient_variation(fignum,fsize,numbins,quantiles)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'drag_coefficient_components_wind_speed.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'drag_coefficient_components_wind_speed.png'],'-dpng',dpi_string);
end

%% Figure 9 - impact on stress estimation via C_D

fignum = 9;
figure(fignum)
set(gcf,'Position',boxpos.*[1 1 2 1])

wind_stress_overestimation(fignum,fsize,numbins,quantiles)

switch print_option
    case 'none'
    case 'svg'
        figure(fignum);print([out_path 'tau_percent_overestimation.svg'],'-dsvg')
    case 'png'
        figure(fignum);print([out_path 'tau_percent_overestimation.png'],'-dpng',dpi_string);
end
