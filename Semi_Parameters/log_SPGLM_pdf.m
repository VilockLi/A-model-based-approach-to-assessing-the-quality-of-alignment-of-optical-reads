function lik=log_SPGLM_pdf(seq,observed_seq,ref_Y,ref_phi,ref_beta,Delta_vec)
% This function calculate the log likelihood of fitted semi-parametric
% distribution

% seq: the original sequence 
% observed_seq: the observed sequence
% ref_Y: 'y' values for fitted semi-parametric distribution 
% ref_phi: fitted phi function (a vector)
% ref_beta: fited beta
% Delta_vec: differences vector of ref_Y

X = sqrt(seq);
X = X(:);
observed_seq = observed_seq(:);
Y = (sqrt(observed_seq)-X) ./ X;
X = [X, ones(length(X), 1)];
phi_vec = [];
for i = 1:length(Y)
    phi_vec = [phi_vec, linear_phi_inter(ref_Y, ref_phi,Y(i))];
end
if sum(phi_vec == -inf) > 0
    lik = -inf;
else
    theta_vec = X * ref_beta;
    n = length(Y);
    term_1 = exp(theta_vec * ref_Y');
    term_2 = exp(repmat(ref_phi', n, 1));
    term_33 = term_1 .* term_2;
    term_3 = repmat(Delta_vec',n,1) ./ (repmat((ref_phi(2:end) - ref_phi(1:end-1))', n, 1) + theta_vec * Delta_vec');
    c = term_33(:, 2:end) - term_33(:, 1:end - 1);
    J = log(sum(c .* term_3, 2));
    lik = sum(phi_vec) + sum(theta_vec .* Y) - sum(J);
end