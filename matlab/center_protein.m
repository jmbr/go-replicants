function Q = center_protein(P)

num_atoms = size(P, 1);
center_of_mass = mean(P, 1);
Q = P - repmat(center_of_mass, num_atoms, 1);
