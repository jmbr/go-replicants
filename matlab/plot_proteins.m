clear all;
close all;

addpath('/home/jmbr/sources/matlab-utils');

P = center_protein(protein_1pgb);
Q = center_protein(protein_2gb1);

[d1 s1] = mean_distance_between_CAs(P);
[d2 s2] = mean_distance_between_CAs(Q);
disp(sprintf('Average distance between alpha carbons in 1PGB: %g (+/- %g) Angstroms.', d1, s1));
disp(sprintf('Average distance between alpha carbons in 2GB1: %g (+/- %g) Angstroms.', d2, s2));

[A b] = compute_affinity(P', Q');
P = transpose(A*transpose(P));

dmax = 10;

% dnat1 = contact_map(P, dmax);
% figure;
% surface(dnat1);
% title('1PGB');
% axis ij;
figure;
dnat2 = contact_map(Q, dmax);
surface(dnat2);
title('2GB1');
axis ij;

a = 0.5;                                % Less than 1 Angstrom.

% simulation.protein = P;
% simulation.num_atoms = size(simulation.protein, 1);
% simulation.dnat = dnat1;
% simulation.a = a;
% simulation.dmax = dmax;

U1 = potential(simulation);
% U2 = potential(P, dnat2, a);

% figure;
% hold on;
% plot3(P(:, 1), P(:, 2), P(:, 3), 'o-b');
% plot3(Q(:, 1), Q(:, 2), Q(:, 3), 'o-r');
% axis off;
% axis equal;
% hold off;
