function lik_hood=log_concave_lik(phi_vec,Y,Delta_vec,theta_vec,term_1,term_2)
% This function returns the likelihood of semi-parametric distribution.

% phi_vec is the non-paramatric part, which is concave in shape
% beta_coe is the coefficient of X
% X is the design matrix
% Y is the response vector
% term_1 and term_2 are terms that will not change during the iteration
N=length(Y);
unique_Y = unique(Y);
% theta_vec is the vector of theta_i, thus has length N
% note theta_vec essentially is the vector of ``predicted values'' of y
% delta_vec is the vector of deltas, thus has length K-1
% The log likelihood
term_33 = term_1 .* term_2;
term_3 = repmat(Delta_vec', N, 1) ./ (repmat((phi_vec(2:end) - phi_vec(1:end-1))', N, 1) + theta_vec * Delta_vec');
c = term_33(:, 2:end) - term_33(:, 1:end - 1);
J = log(sum(c .* term_3, 2));
L = zeros(N,1);
for i = 1:N
    ind = find(unique_Y == Y(i));
    L(i) = L(i) + phi_vec(ind);
end
lik_hood = sum(L) + sum(theta_vec .* Y) - sum(J);
