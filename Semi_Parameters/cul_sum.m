function cul_vec = cul_sum(v)
% this function calculates the culmulative sum

% v: vector of numbers.
cul_vec = [];
sum = 0;
for i = 1:length(v)
    sum = sum + v(i);
    cul_vec = [cul_vec, sum];
end
