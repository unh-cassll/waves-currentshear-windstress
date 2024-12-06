% Given inputs of
% * independent variable in_x
% * dependent variable in_y
% * number of bins (with equal size)
% * quantiles to compute
%
% ... evaluates the quantiles of the dependent variable in each bin
%
% N. Laxague 2024
%
function [x_quantiles,y_quantiles,binsize] = compute_quantiles_fixed_binsize(in_x,in_y,numbins,quantiles)

inds_nan = isnan(in_x) | isnan(in_y);
in_x(inds_nan) = [];
in_y(inds_nan) = [];

binsize = floor(length(in_x)/numbins);

[sorted_x,order] = sort(in_x);
sorted_y = in_y(order);

x_quantile_block = NaN*ones(numbins,length(quantiles));
y_quantile_block = x_quantile_block;

for i = 1:numbins

    ind_s = (i-1)*binsize+1;
    ind_f = i*binsize;
    snipx = sorted_x(ind_s:ind_f);
    sort_snipx = sorted_x(ind_s:ind_f);
    sort_snipy = sort(sorted_y(ind_s:ind_f));

    for j = 1:length(quantiles)
        ind_quantile = floor(quantiles(j)*length(snipx)/100);
        x_quantile_block(i,j) = sort_snipx(ind_quantile);
        y_quantile_block(i,j) = sort_snipy(ind_quantile);
    end

end

x_quantiles = x_quantile_block;
y_quantiles = y_quantile_block;

