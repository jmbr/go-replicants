#!/usr/bin/awk -f

BEGIN {
  num_alpha_carbons = 0;
  printf("const vector positions[] = {");
}

/^ATOM.*CA/ {
  printf("{ %g, %g, %g }, ", $7, $8, $9);
  ++num_alpha_carbons;
}

END {
  printf("};\nconst size_t num_atoms = %d;  // Number of alpha carbons.\n", num_alpha_carbons);
}
