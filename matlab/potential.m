function U = potential(simulation)

if simulation.a <= 0
    error('Invalid value for a. It should be strictly positive.');
end

U = -simulation.num_atoms;

for i = 1:simulation.num_atoms
    for j = i+1:simulation.num_atoms
        U = U + pairwise_potential(simulation, i, j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = pairwise_potential(simulation, i, j)

w_ij = 1/simulation.a^2;

r = distance_between_amino_acids(simulation.protein, i, j);

if abs(r - simulation.dnat(i, j)) < simulation.a
    u = w_ij*((r - simulation.dnat(i, j))^2 - simulation.a^2);
else
    u = 0;
end
