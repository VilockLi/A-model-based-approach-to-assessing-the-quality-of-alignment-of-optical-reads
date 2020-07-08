function [log_post_lik, post_part, prior_part, candi_allignment]=appro_log_post_vec(ori_b,ori_e,ori_seq,lam,L,observed_seq,T,ref_Y,ref_phi,ref_beta)
% this function calculate the approximation of log posterior likelihood of 
% the semi-parametric distribution

% b,e: start and end position of original sequence 
% seq: original sequence 
% lam: lambda for possion process (cutting procedure) 
% L: physical length of whole genome 

c = (L - ori_e + ori_b + 2) / L;
% We take log scale up to a constant
log_prior = lam * c;
Delta_vec = ref_Y(2:end) - ref_Y(1:end - 1);
candi_allignment = [];
likelihood_allignment = [];
observed_seq = observed_seq(:);
ori_seq = ori_seq(:);
relMat = abs((repmat(sqrt(observed_seq), 1, length(ori_seq)) - repmat(sqrt(ori_seq), 1, length(observed_seq))')...
    ./repmat(sqrt(ori_seq), 1, length(observed_seq))');
scoreMat = inf * ones(size(relMat));
recordCell = cell(length(observed_seq), length(ori_seq));
for i = 1: length(ori_seq) - length(observed_seq) + 1
    scoreMat(1, i) = relMat(1, i);
    recordCell{1, i} = [i];
end
for i = 2:length(observed_seq)
    for j = i: i + length(ori_seq) - length(observed_seq)
        ind = find(scoreMat(i-1, 1:j-1) ==  min(scoreMat(i-1, 1:j-1)));
        ind = ind(1);
        scoreMat(i, j) = relMat(i, j) + scoreMat(i-1, ind);
        recordCell{i, j} = [recordCell{i-1, ind}, j];
    end
end
allignInd = recordCell{length(observed_seq), scoreMat(end, :) == min(scoreMat(end,:))};

lik_temp = exp(log_SPGLM_pdf(ori_seq(allignInd)',observed_seq,ref_Y,ref_phi,ref_beta,Delta_vec));
adjInd = setdiff(1:length(ori_seq), allignInd);
adjInd(adjInd > allignInd(end)) = [];
adjInd(adjInd < allignInd(1)) = [];
for i = 1:length(adjInd)
    lik_temp = lik_temp * concave_den_cdf(ref_Y,ref_phi,ref_beta,Delta_vec,T, ori_seq(adjInd(i)));
end


log_post_lik = lam*c + log(lik_temp) - lam+2*log(lam/L);
post_part = log(lik_temp);
prior_part = lam * c;