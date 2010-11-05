#!/usr/bin/awk -f

BEGIN {
  num_alpha_carbons = 0;
  printf("P = [");
}

/^ATOM.*CA/ {
  printf("%3.3g, %3.3g, %3.3g;", $7, $8, $9);
  ++num_alpha_carbons;
}

END {
  printf("];\n");
}
