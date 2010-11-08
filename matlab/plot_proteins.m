clear all;
close all;

addpath('/home/jmbr/sources/matlab-utils');

P = center_protein(protein_1pgb);
Q = center_protein(protein_2gb1);

[d1 s1] = mean_distance_between_CAs(P);
[d2 s2] = mean_distance_between_CAs(Q);
disp(sprintf('Average distance between alpha carbons: %g (+/- %g) Angstroms.', d1, s1));
disp(sprintf('Average distance between alpha carbons: %g (+/- %g) Angstroms.', d2, s2));

[A b] = compute_affinity(P', Q');
P = transpose(A*transpose(P));

hold on;
plot3(P(:, 1), P(:, 2), P(:, 3), 'o-b');
plot3(Q(:, 1), Q(:, 2), Q(:, 3), 'o-r');
axis off;
axis equal;
hold off;
