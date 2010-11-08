function dnat = contact_map(P, dmax)

if nargin == 1
    dmax = +Inf;
end

N = size(P, 1);
dnat = zeros(N);
for i = 1:N
    for j = i+1:N
        d = distance_between_amino_acids(P, i, j);
        if abs(d) <= dmax
            dnat(i, j) = d;
        else
            dnat(i, j) = +Inf;
        end
        dnat(j, i) = dnat(i, j);
    end
end
