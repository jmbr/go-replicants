function [d s] = mean_distance_between_CAs(P)

num_atoms = size(P, 1);
Q = zeros(1, num_atoms-1);
for k = 1:num_atoms-1
    Q(k) = norm(P(k+1, :) - P(k, :));
end

d = mean(Q);
s = std(Q);
