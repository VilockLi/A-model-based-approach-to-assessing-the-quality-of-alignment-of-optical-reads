function A_val=cal_A_val(ref_Y,ref_beta,ref_phi,Delta_vec,x)
% this function calculates the A value (normalization part) of our
% distribtuion

% ref_Y, ref_beta, ref_phi: all fitted from semi-parametric distribution
% Delta_vec: differences vector of ref_Y
% x: corresponding x

theta_vec = x * ref_beta;
n = length(x(:, 1));
term_1 = exp(theta_vec * ref_Y');
term_2 = exp(repmat(ref_phi', n, 1));
term_33 = term_1 .* term_2;
term_3 = repmat(Delta_vec',n,1) ./ (repmat((ref_phi(2:end) - ref_phi(1:end - 1))', n, 1) + theta_vec * Delta_vec');
c = term_33(:,2:end) - term_33(:, 1:end - 1);
A_val = log(sum(c .* term_3, 2));
