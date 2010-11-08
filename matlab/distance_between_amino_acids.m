function r = distance_between_amino_acids(protein, i, j)

triple_scalar_product = @(u, v, w) dot(u, cross(v, w));

if abs(i - j) == 3
    u = protein(i+1, :) - protein(i, :);
    v = protein(i+2, :) - protein(i+1, :);
    w = protein(i+3, :) - protein(i+2, :);
    signum = sign(triple_scalar_product(protein(i, :), protein(i+1, :), protein(i+2, :)));
else
    signum = 1;
end
        
r = signum*norm(protein(i, :) - protein(j, :));
