addpath Semi_Parameters
load phi_opt
load blocks
load opt_beta

%% fitting linear spline to phi
% to speed up the calculation, we only keep few points to approximate the
% phi() function
slope = [];
sort_Y = sort(Y);
% calculate the slopes
for i = 2:length(Y)
    slope = [slope,(phi_opt(i) - phi_opt(i - 1)) / (sort_Y(i) - sort_Y(i - 1))];
end
% calculate the relative difference between slopes
rel_dif = [];
for i = 2:length(slope)
    rel_dif = [rel_dif,abs((slope(i) - slope(i - 1))/slope(i - 1))];
end
% set the threshold, we only consider points that lead to a change of
% relative difference larger than the threshold
tre_slope = 0.02;
fit_ind = [];
for i = find(rel_dif > tre_slope)
    fit_ind = [fit_ind, i, i + 1];
end
fit_ind = [1, fit_ind, length(sort_Y)];
fit_ind = unique(fit_ind);
% randomly pick several points
v_1 = fit_ind(fit_ind > 1500 & fit_ind < 2000);
v_2 = fit_ind(fit_ind > 2000);
fit_ind = floor(unique([1,quantile(v_1, 0.25), quantile(v_1, 0.5),...
    quantile(v_1,0.75), quantile(v_2, 0.2)...
    quantile(v_2,0.4), quantile(v_2, 0.6), quantile(v_2, 0.8), quantile(v_2, 0.9),...
    length(sort_Y)]));
ref_Y = sort_Y(fit_ind);
ref_phi = phi_opt(fit_ind);

cd Semi_Parameters
save('approximated_phi', 'ref_phi')
save('approximated_y', 'ref_Y')
cd ..