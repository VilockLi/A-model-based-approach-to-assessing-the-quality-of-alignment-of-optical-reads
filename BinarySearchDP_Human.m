clear all; clc
addpath Synthetic_Data/HumanReference
addpath Semi_Parameters

%% load optical mapping data file
genome_id = 1;
outputFilename = 'output';
load(strcat('fake_genome_', num2str(genome_id)))
load(strcat('fake_reads_', num2str(genome_id)))

%% load algorithm parameter file
load approximated_phi
load approximated_y
load opt_beta


%% pre-processing
allReads = fake_reads_1;
numReads = size(allReads, 1);
ori_ref = refGenome';
%% parameters
aveLength = 250000; T = 400; L = sum(ori_ref); lam = floor(L / aveLength);
cnt_bnd = 1.5;

ref_beta = beta_opt;
% quantities used in calculating likelihood
Delta_vec = ref_Y(2:end) - ref_Y(1:end - 1);
delta_vec = ref_phi(2:end) - ref_phi(1:end - 1);


% anything below T is not observable
% find the cutoff points
for x = 10:10:2*T
    prob = concave_den_cdf(ref_Y,ref_phi,ref_beta,Delta_vec,T,x);
    if prob < 0.99
        break
    end
end
% any fragment shorter than cutoff_low is likely to be unobservable
cutoff_low = x - 10;
% any fragment longer than cutoff_up is very likely to be observable
cutoff_up = (sqrt(T) / (1 - abs(min(ref_Y))))^2;
% support
support_ref_Y = [min(ref_Y), max(ref_Y)];


for id_num_loop = 1
    % the id of the read
    id_num = id_num_loop;
    observed_seq = allReads(id_num, :);
    % delete the zeros entries
    observed_seq(1:2) = [];
    observed_seq(observed_seq == 0) = [];
    observed_seq(isnan(observed_seq)) = [];
    T_observed_seq = observed_seq;
    m = length(T_observed_seq);
    % remove invalid reads
    if m < 4
        all_data = [nan];
        cd(strcat('oldResultsSEMI/genome_', num2str(genome_id)))
        file_name = strcat(outputFilename, num2str(id_num));
        save(file_name, 'all_data');
        cd ../..
        continue
    end
    sqrt_observed_seq = sqrt(T_observed_seq);
    
    % the threshold we use to conduct prefiltering
    cnt_tred = 0.375;
    % we conduct the filtering through several rounds
    cell_pos_temp  = [];
    all_data = [];
    % the stopping criteria
    while cnt_tred <= cnt_bnd - 1.125
        cell_pos = [];
        ref = ori_ref;
        % the filtered genome
        ref(ref < T*cnt_tred) = [];
        ref_filter = ref(ref > T*cnt_tred);
        % the positions of fragments that cannot be dropped
        cant_drop_ind = ref_filter > cutoff_up;
        
        Bool_mat = [];
        for i = 1:length(T_observed_seq)
            vec = ((sqrt_observed_seq(i) - sqrt(ref_filter)) ./ sqrt(ref_filter))';
            if i == 1
                vec(end - length(T_observed_seq) + 1 + i:end) = [];
            elseif i < length(T_observed_seq)
                vec(1:i - 1) = [];
                vec(end - length(T_observed_seq) + 1 + i:end) = [];
            else
                vec(1:i - 1) = [];
            end
            vec = vec(:)';
            Bool_mat = [Bool_mat; support_ref_Y(1) <= vec & vec <= support_ref_Y(2)];
        end
        % first round: k = m
        ind_1 = find(Bool_mat(1, :)~=0);
        check_1 = sum(Bool_mat(:, ind_1));
        feasible_1 = ind_1(check_1 == length(T_observed_seq));
        cell_pos = [feasible_1', feasible_1' + length(T_observed_seq) - 1];
        
        % second round: k = m + 1
        % drop can happen at 1 out of (m-1) positions, we do not allow it to drop
        % at the end of sequence (which basically is the first case)
        mat_1 = nchoosek(1:length(T_observed_seq) - 1, 1);
        survive_ind = intersect(ind_1, 1:length(Bool_mat) - 1);
        for i = 1 : size(mat_1, 1)
            % drop happens between ith and (i + 1)th position
            temp_mat = [Bool_mat(1:i, survive_ind); Bool_mat(i + 1:end, 1 + survive_ind)];
            temp_cant_drop = cant_drop_ind(survive_ind + i);
            if i == 1
                % feasible starting points at this round
                cnt = (all(temp_mat)' == 1)&(temp_cant_drop ~= 1);
            else
                % unfeasible starting points at this round
                previous_check = sum(temp_mat(1:i,:)) >= i;
                del_ind = find(previous_check'+cnt == 0);
                % update survive starting points
                survive_ind(del_ind) = [];
                cnt(del_ind) = [];
                temp_mat(:, del_ind) = [];
                temp_cant_drop(del_ind) = [];
                cnt = cnt + ((all(temp_mat)' == 1) & (temp_cant_drop == 0));
            end
            ind_temp = ((all(temp_mat)' == 1) & (temp_cant_drop == 0));
        end
        feasible_2 = survive_ind(cnt ~= 0);
        cell_pos = [cell_pos; feasible_2', feasible_2' + length(T_observed_seq)];
        display('finish mat_1')
        
        % third round: k = m + 2
        % drop can happen at 2 out of m positions, we do not allow it to drop
        % at the end of sequence (which basically is the first  or second case)
        mat_2 = nchoosek(1:length(T_observed_seq), 2);
        dif_2 = mat_2(:,2) - mat_2(:, 1);
        survive_ind = intersect(ind_1, 1:length(Bool_mat) - 2);
        for i = 1 : size(mat_2,1)
            if mat_2(i,2) - mat_2(i,1) == 1
                % this means we dropped two adjacent fragments
                temp_mat = [Bool_mat(1:mat_2(i,1), survive_ind);...
                    Bool_mat(mat_2(i,1) + 1:end, 2 + survive_ind)];
                temp_cant_drop = cant_drop_ind(survive_ind+mat_2(i, 1))+...
                    cant_drop_ind(survive_ind + mat_2(i,1) + 1);
            else
                % we have gap between two dropped fragments
                temp_mat = [Bool_mat(1:mat_2(i, 1), survive_ind);...
                    Bool_mat(mat_2(i, 1) + 1:mat_2(i, 2) - 1, 1 + survive_ind);...
                    Bool_mat(mat_2(i, 2):end, 2 + survive_ind)];
                temp_cant_drop = cant_drop_ind(survive_ind + mat_2(i, 1)) +...
                    cant_drop_ind(survive_ind + mat_2(i, 2));
            end
            if i == 1
                % feasible starting points at this round
                cnt = (all(temp_mat)' == 1)&(temp_cant_drop == 0);
            else
                % unfeasible starting points at this round
                previous_check = sum(temp_mat(1:mat_2(i, 1), :)) >= mat_2(i, 1);
                del_ind = find(previous_check' + cnt == 0);
                % update survive starting points
                survive_ind(del_ind) = [];
                cnt(del_ind) = [];
                temp_mat(:, del_ind) = [];
                temp_cant_drop(del_ind) = [];
                cnt = cnt + ((all(temp_mat)' == 1) & (temp_cant_drop == 0));
            end
            ind_temp = ((all(temp_mat)' == 1) & (temp_cant_drop == 0));
        end
        feasible_3 = survive_ind(cnt ~= 0);
        cell_pos = [cell_pos;feasible_3', feasible_3' + length(T_observed_seq) + 1];
        display('finish mat_2')
        
        % the fourth round k = m+3
        mat_3 = nchoosek(1:length(T_observed_seq) + 1, 3);
        dif_3 = mat_3(:, 3) - mat_3(:, 1);
        survive_ind = intersect(ind_1, 1:length(Bool_mat) - 3);
        for i = 1 : size(mat_3, 1)
            if mat_3(i, 2)-mat_3(i, 1) == 1 & mat_3(i, 3) - mat_3(i, 2) == 1
                % this means we dropped two adjacent fragments
                temp_mat = [Bool_mat(1:mat_3(i, 1), survive_ind);...
                    Bool_mat(mat_3(i, 1) + 1:end,3 + survive_ind)];
                temp_cant_drop = cant_drop_ind(survive_ind + mat_3(i, 1)) +...
                    cant_drop_ind(survive_ind+mat_3(i, 1)+1) +...
                    cant_drop_ind(survive_ind+mat_3(i, 1)+2);
            elseif mat_3(i, 2) - mat_3(i, 1) == 1
                % we have gap between two dropped fragments
                temp_mat = [Bool_mat(1:mat_3(i, 1), survive_ind);...
                    Bool_mat(mat_3(i, 1) + 1:mat_3(i, 3) - 2, 2 + survive_ind);...
                    Bool_mat(mat_3(i, 3) - 1:end, 3 + survive_ind)];
                temp_cant_drop = cant_drop_ind(survive_ind + mat_3(i, 1)) +...
                    cant_drop_ind(survive_ind + mat_3(i, 1) + 1) +...
                    cant_drop_ind(survive_ind + mat_3(i, 3));
            elseif mat_3(i,3) - mat_3(i,2) == 1
                temp_mat = [Bool_mat(1:mat_3(i, 1), survive_ind);...
                    Bool_mat(mat_3(i, 1)+1:mat_3(i, 2) - 1, 1 + survive_ind);...
                    Bool_mat(mat_3(i, 2):end, 3 + survive_ind)];
                temp_cant_drop = cant_drop_ind(survive_ind+mat_3(i, 1)) +...
                    cant_drop_ind(survive_ind + mat_3(i, 2)) +...
                    cant_drop_ind(survive_ind + mat_3(i,3));
            else
                temp_mat = [Bool_mat(1:mat_3(i, 1), survive_ind);...
                    Bool_mat(mat_3(i, 1)+1:mat_3(i,2) - 1, 1 + survive_ind);...
                    Bool_mat(mat_3(i, 2):mat_3(i,3) - 2, 2 + survive_ind);...
                    Bool_mat(mat_3(i, 3) - 1:end, 3 + survive_ind)];
                temp_cant_drop=cant_drop_ind(survive_ind+mat_3(i,1))+...
                    cant_drop_ind(survive_ind + mat_3(i, 2)) +...
                    cant_drop_ind(survive_ind + mat_3(i, 3));
            end
            if i == 1
                % feasible starting points at this round
                cnt = (all(temp_mat)' == 1) & (temp_cant_drop == 0);
            else
                % unfeasible starting points at this round
                previous_check = sum(temp_mat(1:mat_3(i, 1), :)) >= mat_3(i, 1);
                del_ind = find(previous_check' + cnt == 0);
                % update survive starting points
                survive_ind(del_ind) = [];
                cnt(del_ind) = [];
                temp_mat(:, del_ind) = [];
                temp_cant_drop(del_ind) = [];
                cnt = cnt + ((all(temp_mat)' == 1) & (temp_cant_drop == 0));
            end
            ind_temp = ((all(temp_mat)' == 1) & (temp_cant_drop == 0));
        end
        feasible_4 = survive_ind(cnt ~= 0);
        cell_pos = [cell_pos; feasible_4', feasible_4' + length(T_observed_seq) + 2];
        display('finish mat_3')
        if ~isempty(cell_pos)
            new_id = find(ori_ref > cnt_tred*T);
            cell_pos (cell_pos(:,1)==1,:) = [];
            l_of_cell = size(cell_pos, 1);
            cell_pos = [new_id(cell_pos(:, 1)),new_id(cell_pos(:, 2))];
        end
        cnt_tred = cnt_tred + 0.375;
        cell_pos_temp = [cell_pos_temp; cell_pos];
    end
    %% Calculate the likelihoods
    ref = ori_ref;
    cul_ref = cul_sum(ref);
    cell_pos = unique(cell_pos_temp, 'rows');
    if ~isempty(cell_pos)
        [l_of_cell, w_of_cell] = size(cell_pos);
        if w_of_cell == 1
            cell_pos = cell_pos';
            l_of_cell = 1;
        end
        prob = zeros(1,l_of_cell);
        for i = 1:l_of_cell
            ori_e = cul_ref(cell_pos(i,2));
            % get the sequence
            ori_seq = ref(cell_pos(i,1):cell_pos(i,2));
            ori_b = cul_ref(cell_pos(i,1) - 1);
            post = appro_log_post_vec(ori_b, ori_e, ori_seq, lam, L, observed_seq, T, ref_Y,...
                ref_phi, ref_beta);
            prob(i) = post;
        end
        % remove the -inf ones
        prob_temp = prob;
        prob_temp(prob_temp == -inf) = [];
        if ~isempty(prob_temp)
            prob(prob == -inf) = min(prob_temp) - 100;
        else
            prob(prob == -inf) = - 100;
        end
        [~, ind] = sort(prob, 'descend');
        un_norm = sum(exp(prob));
        prob_1 = exp(prob(ind(1))) / un_norm;
        start_1 = cell_pos(ind(1), 1);
        end_1 = cell_pos(ind(1), 2);
        all_data = zeros(length(prob), 7);
        all_data(:, 5:6) = cell_pos(ind, :);
        all_data(:, 7) = exp(prob(ind)) / un_norm;
        all_data(1, 1) = id_num;
        all_data(1, 2) = un_norm;
        all_data(1, 3) = length(observed_seq);
        all_data(1, 4) = l_of_cell;
    end
    id_num_loop
end

%% The candidates are stored in all_data, ranked by their corresponding posterior likelihoods.
all_data = array2table(all_data);
all_data.Properties.VariableNames(1:end) = {'id', 'normalizing_constant',...
    'length_of_read', 'number_of_candidates', 'start_pos', 'end_pos', 'posterior_likelihood'};