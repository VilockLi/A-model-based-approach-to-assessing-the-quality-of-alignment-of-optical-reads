function phi_val=linear_phi_inter(ref_Y,ref_phi,y)
% this function finds the linear interpolation of y, i.e. \phi(y)

% ref_Y: 'y' values for fitted semi-parametric distribution 
% ref_phi: fitted phi function (a vector)
% y: the point we want to evaluate on


if y < ref_Y(1) || y >  ref_Y(end)
    phi_val = -inf;
else
    ind = find(ref_Y == y);
    if isempty(ind)
        ind_1 = sum(ref_Y < y);
        ind_2 = ind_1 + 1;
        phi_val = (ref_phi(ind_2) - ref_phi(ind_1)) / (ref_Y(ind_2) - ref_Y(ind_1)) * y +...
            (ref_phi(ind_1) * ref_Y(ind_2) - ref_phi(ind_2) * ref_Y(ind_1)) / (ref_Y(ind_2) - ref_Y(ind_1));
    else
        phi_val = ref_phi(ind);
    end
end
