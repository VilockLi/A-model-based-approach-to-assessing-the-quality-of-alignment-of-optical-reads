function [gen_seq, ori_seq]=gen_seq_semi(seq, ref_Y, ref_beta, ref_phi, Delta_vec, T)
% this function generates mapping based on semi-parametric distribution
seq = seq(:);
m = length(seq);
rand_vec = zeros(m,1);
for i = 1:m
    x = [sqrt(seq(i)),1];
    A_val = cal_A_val(ref_Y,ref_beta,ref_phi,Delta_vec,x);
    theta_i = sum(x'.*ref_beta);
    c = -inf;
    for j = 1:length(ref_Y)
        y = ref_Y(j);
        phi_val = linear_phi_inter(ref_Y,ref_phi,y);
        f_val = exp(theta_i*y+phi_val-A_val);
        if f_val > c
            c = f_val;
            best_y = y;
        end
    end
    a = min(ref_Y);
    b = max(ref_Y);
    accept = 0;
    while ~accept
        U_1 = rand(1);
        U_2 = rand(1);
        r_1 = U_1 *(b-a)+a;
        r_2 = U_2 *c;
        phi_val = linear_phi_inter(ref_Y,ref_phi,r_1);
        f_val = exp(theta_i*r_1+phi_val-A_val);
        if f_val > r_2
            accept = 1;
            rand_vec(i) = r_1;
        end
    end
end   
gen_seq = sqrt(seq).*(ones(m,1)+rand_vec);
gen_seq = gen_seq.^2;
idx = find(gen_seq < T);
seq(idx) = nan;
ori_seq = seq;
gen_seq(idx) = nan;
