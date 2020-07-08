clc; clear all 
%% this script generates new synthetic mappings 
addpath '../Synthetic_Data'
addpath '../Semi_Parameters'
addpath '../Fit_Semi'

load('fragment_len_1')
load('fragment_len_2')
load('fragment_len_3')
load('dictionary_genome')
load('phyLength')
load('len_vec')
load('opt_phi')
load('blocks')
load('opt_beta')


%% parameters
% to sample from semi-parametric distribution, we do not use the
% approximated distribution
Y = sort(Y);
Delta_vec = Y(2:end)-Y(1:end-1);
delta_vec = phi_opt(2:end)-phi_opt(1:end-1);

support_Y = [min(Y),max(Y)];

% the threshold
T = 400;
% the number of sequences
n = 1000;

mu = 0;
sigma = 0.1;
semi_Y = Y;

% we only consider the first 23 genomes 
ind_1 = find(frag_len_1 == 0);
refGenome = frag_len_1(1:ind_1(24));
refGenome(ind_1(1:24)) = [];
L = sum(refGenome);
aveLength = mean(phyLength);
cumref = cumsum(refGenome);

% now we shuffle the genome
% len_vec_1 = randsample(len_vec,n);
% start_point_1 = randsample(length(refGenome)-max(len_vec-1),n);
% end_point_1 = start_point_1 + len_vec_1-1;

% the first two columns are start/end positions. 
fake_original = zeros(n, 2 + max(len_vec));
for i = 1:n
    segStart = 0;
    segStop = 0;
    while segStart + 4 >= segStop | segStart + max(len_vec) <= segStop
        queryLength = phyLength(randsample(n, 1));
        phyStart = randsample(floor(L / 1000), 1) * 100 + randsample(1000, 1);
        phyStop = min(floor(phyStart + exprnd(aveLength)), L);
        if sum(phyStart == cumref)
            phyStart = phyStart - 1;
        end
        if sum(phyStop == cumref)
            phyStop = phyStop + 1;
        end
        segStart = sum(cumref < phyStart);
        segStop = sum(cumref < phyStop);
    end
    fake_original(i,1) = segStart;
    fake_original(i,2) = segStop;
    fake_original(i, 3:segStop - segStart + 3) = refGenome(segStart:segStop);
end

fake_read_normal_1 = zeros(n, max(len_vec));
for i = 1:n
    new_vec = gen_seq_norm(refGenome(fake_original(i,1):fake_original(i,2)), mu, sigma, T);
    fake_read_normal_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = floor(new_vec);
end

fake_read_semi_1 =  zeros(n, max(len_vec));
for i = 1:n
    seq = refGenome(fake_original(i,1):fake_original(i,2));
    new_vec = gen_seq_semi(seq, Y, beta_opt, phi_opt, Delta_vec, T);
    fake_read_semi_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = floor(new_vec);
end

fake_read_nb_1 = zeros(n, max(len_vec));
for i = 1:n
    new_vec = nbinrnd(refGenome(fake_original(i,1):fake_original(i,2)), 0.5);
    fake_read_nb_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = new_vec;
end

% put all these reads together and delete the bad ones 
fake_reads_1 = zeros(3 * n, 2 + max(len_vec));
fake_reads_1(:, 3:end) = [fake_read_normal_1; fake_read_semi_1; fake_read_nb_1];
fake_reads_1(1:n, 1) = 1:n;
fake_reads_1(n + 1:2 * n, 1) = 1:n;
fake_reads_1(2 * n + 1: 3 * n, 1) = 1:n;
fake_reads_1(1:n, 2) = 1;
fake_reads_1(n + 1:2 * n, 2) = 2;
fake_reads_1(2 * n + 1:3 * n, 2) = 3;
del_ind = [];
for i = 1:3*n
    temp_vec = fake_reads_1(i,3:end);
    temp_vec(temp_vec<T) = [];
    if length(temp_vec)>0
        fake_reads_1(i,3:(2+length(temp_vec))) = temp_vec;
        fake_reads_1(i,(3+length(temp_vec)):end) = 0;
    end
    if sum(temp_vec < 0) > 0 
        del_ind = [del_ind,i];
    end
end
if ~isempty(del_ind)
    fake_reads_1(del_ind,:) = [];
end

cd '../Synthetic_Data'
save('fake_genome_1','refGenome')
save('fake_original_1','fake_original')
save('fake_read_normal_1','fake_read_normal_1')
save('fake_read_semi_1','fake_read_semi_1')
save('fake_read_nb_1','fake_read_nb_1')
save('fake_reads_1','fake_reads_1')
cd ../Generate_Synthetic


%% genome 2
% we only consider the first 23 genomes 
% we only consider the first 23 genomes 
ind_1 = find(frag_len_2 == 0);
refGenome = frag_len_2(1:ind_1(24));
refGenome(ind_1(1:24)) = [];
L = sum(refGenome);
aveLength = mean(phyLength);
cumref = cumsum(refGenome);

% now we shuffle the genome
% len_vec_1 = randsample(len_vec,n);
% start_point_1 = randsample(length(refGenome)-max(len_vec-1),n);
% end_point_1 = start_point_1 + len_vec_1-1;

% the first two columns are start/end positions. 
fake_original = zeros(n, 2 + max(len_vec));
for i = 1:n
    segStart = 0;
    segStop = 0;
    while segStart + 4 >= segStop | segStart + max(len_vec) <= segStop
        queryLength = phyLength(randsample(n, 1));
        phyStart = randsample(floor(L / 1000), 1) * 100 + randsample(1000, 1);
        phyStop = min(floor(phyStart + exprnd(aveLength)), L);
        if sum(phyStart == cumref)
            phyStart = phyStart - 1;
        end
        if sum(phyStop == cumref)
            phyStop = phyStop + 1;
        end
        segStart = sum(cumref < phyStart);
        segStop = sum(cumref < phyStop);
    end
    fake_original(i,1) = segStart;
    fake_original(i,2) = segStop;
    fake_original(i, 3:segStop - segStart + 3) = refGenome(segStart:segStop);
end

fake_read_normal_1 = zeros(n, max(len_vec));
for i = 1:n
    new_vec = gen_seq_norm(refGenome(fake_original(i,1):fake_original(i,2)), mu, sigma, T);
    fake_read_normal_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = floor(new_vec);
end

fake_read_semi_1 =  zeros(n, max(len_vec));
for i = 1:n
    seq = refGenome(fake_original(i,1):fake_original(i,2));
    new_vec = gen_seq_semi(seq, Y, beta_opt, phi_opt, Delta_vec, T);
    fake_read_semi_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = floor(new_vec);
end

fake_read_nb_1 = zeros(n, max(len_vec));
for i = 1:n
    new_vec = nbinrnd(refGenome(fake_original(i,1):fake_original(i,2)), 0.5);
    fake_read_nb_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = new_vec;
end

% put all these reads together and delete the bad ones 
fake_reads_1 = zeros(3 * n, 2 + max(len_vec));
fake_reads_1(:, 3:end) = [fake_read_normal_1; fake_read_semi_1; fake_read_nb_1];
fake_reads_1(1:n, 1) = 1:n;
fake_reads_1(n + 1:2 * n, 1) = 1:n;
fake_reads_1(2 * n + 1: 3 * n, 1) = 1:n;
fake_reads_1(1:n, 2) = 1;
fake_reads_1(n + 1:2 * n, 2) = 2;
fake_reads_1(2 * n + 1:3 * n, 2) = 3;
del_ind = [];
for i = 1:3*n
    temp_vec = fake_reads_1(i,3:end);
    temp_vec(temp_vec<T) = [];
    if length(temp_vec)>0
        fake_reads_1(i,3:(2+length(temp_vec))) = temp_vec;
        fake_reads_1(i,(3+length(temp_vec)):end) = 0;
    end
    if sum(temp_vec < 0) > 0 
        del_ind = [del_ind,i];
    end
end
if ~isempty(del_ind)
    fake_reads_1(del_ind,:) = [];
end

cd '../Synthetic_Data'
save('fake_genome_2','refGenome')
save('fake_original_2','fake_original')
save('fake_read_normal_2','fake_read_normal_1')
save('fake_read_semi_2','fake_read_semi_1')
save('fake_read_nb_2','fake_read_nb_1')
save('fake_reads_2','fake_reads_1')
cd ../Generate_Synthetic
%% genome 3
ind_1 = find(frag_len_3 == 0);
refGenome = frag_len_3(1:ind_1(24));
refGenome(ind_1(1:24)) = [];
L = sum(refGenome);
aveLength = mean(phyLength);
cumref = cumsum(refGenome);

% now we shuffle the genome
% len_vec_1 = randsample(len_vec,n);
% start_point_1 = randsample(length(refGenome)-max(len_vec-1),n);
% end_point_1 = start_point_1 + len_vec_1-1;

% the first two columns are start/end positions. 
fake_original = zeros(n, 2 + max(len_vec));
for i = 1:n
    segStart = 0;
    segStop = 0;
    while segStart + 4 >= segStop | segStart + max(len_vec) <= segStop
        queryLength = phyLength(randsample(n, 1));
        phyStart = randsample(floor(L / 1000), 1) * 100 + randsample(1000, 1);
        phyStop = min(floor(phyStart + exprnd(aveLength)), L);
        if sum(phyStart == cumref)
            phyStart = phyStart - 1;
        end
        if sum(phyStop == cumref)
            phyStop = phyStop + 1;
        end
        segStart = sum(cumref < phyStart);
        segStop = sum(cumref < phyStop);
    end
    fake_original(i,1) = segStart;
    fake_original(i,2) = segStop;
    fake_original(i, 3:segStop - segStart + 3) = refGenome(segStart:segStop);
end
fake_read_normal_1 = zeros(n, max(len_vec));
for i = 1:n
    new_vec = gen_seq_norm(refGenome(fake_original(i,1):fake_original(i,2)), mu, sigma, T);
    fake_read_normal_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = floor(new_vec);
end

fake_read_semi_1 =  zeros(n, max(len_vec));
for i = 1:n
    seq = refGenome(fake_original(i,1):fake_original(i,2));
    new_vec = gen_seq_semi(seq, Y, beta_opt, phi_opt, Delta_vec, T);
    fake_read_semi_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = floor(new_vec);
end


fake_read_nb_1 = zeros(n, max(len_vec));
for i = 1:n
    new_vec = nbinrnd(refGenome(fake_original(i,1):fake_original(i,2)), 0.5);
    fake_read_nb_1(i,1:fake_original(i,2) - fake_original(i,1) + 1) = new_vec;
end

% put all these reads together and delete the bad ones 
fake_reads_1 = zeros(3 * n, 2 + max(len_vec));
fake_reads_1(:, 3:end) = [fake_read_normal_1; fake_read_semi_1; fake_read_nb_1];
fake_reads_1(1:n, 1) = 1:n;
fake_reads_1(n + 1:2 * n, 1) = 1:n;
fake_reads_1(2 * n + 1: 3 * n, 1) = 1:n;
fake_reads_1(1:n, 2) = 1;
fake_reads_1(n + 1:2 * n, 2) = 2;
fake_reads_1(2 * n + 1:3 * n, 2) = 3;
del_ind = [];
for i = 1:3*n
    temp_vec = fake_reads_1(i,3:end);
    temp_vec(temp_vec<T) = [];
    if length(temp_vec)>0
        fake_reads_1(i,3:(2+length(temp_vec))) = temp_vec;
        fake_reads_1(i,(3+length(temp_vec)):end) = 0;
    end
    if sum(temp_vec < 0) > 0 
        del_ind = [del_ind,i];
    end
end
if ~isempty(del_ind)
    fake_reads_1(del_ind,:) = [];
end

cd '../Synthetic_Data'
save('fake_genome_3','refGenome')
save('fake_original_3','fake_original')
save('fake_read_normal_3','fake_read_normal_1')
save('fake_read_semi_3','fake_read_semi_1')
save('fake_read_nb_3','fake_read_nb_1')
save('fake_reads_3','fake_reads_1')
cd ../Generate_Synthetic