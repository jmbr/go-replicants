struct protein *new_protein_1pgb(void)
{
        const size_t num_atoms = 56;
        const double positions[][3] = {{ 13.935, 18.529, 29.843 }, { 13.088, 19.661, 26.283 }, { 12.726, 17.033, 23.612 }, { 12.179, 17.659, 19.887 }, { 10.253, 15.79, 17.221 }, { 11.082, 16.103, 13.475 }, { 8.009, 15.163, 11.389 }, { 8.628, 13.975, 7.913 }, { 5.213, 12.642, 6.966 }, { 3.589, 12.601, 3.497 }, { 1.291, 15.471, 4.538 }, { 2.624, 17.021, 7.787 }, { 6.24, 18.276, 7.903 }, { 8.146, 20.352, 10.511 }, { 9.429, 20.303, 14.17 }, { 7.799, 20.583, 17.589 }, { 9.057, 20.314, 21.095 }, { 7.744, 19.307, 24.468 }, { 8.907, 19.331, 28.095 }, { 8.985, 15.93, 29.75 }, { 10.637, 14.408, 32.787 }, { 12.217, 11.676, 30.658 }, { 12.452, 10.118, 27.226 }, { 9.683, 7.537, 27.584 }, { 7.249, 10.206, 28.617 }, { 8.146, 12.442, 25.683 }, { 7.82, 9.45, 23.25 }, { 4.309, 8.916, 24.682 }, { 3.067, 12.375, 23.977 }, { 4.843, 12.494, 20.533 }, { 3.247, 9.234, 19.455 }, { -0.277, 10.493, 20.449 }, { 0.341, 13.738, 18.651 }, { 1.53, 11.949, 15.49 }, { -1.401, 9.62, 15.646 }, { -3.917, 12.441, 16.16 }, { -2.542, 14.061, 13.009 }, { -2.684, 10.787, 11.031 }, { 1.109, 10.312, 10.729 }, { 2.165, 6.704, 10.736 }, { 5.826, 6.025, 10.072 }, { 9.219, 4.653, 11.307 }, { 10.629, 6.26, 14.43 }, { 14.117, 6.981, 15.677 }, { 15.402, 8.494, 18.951 }, { 18.798, 10.067, 19.525 }, { 19.598, 10.679, 23.19 }, { 22.584, 12.858, 22.345 }, { 20.338, 15.595, 20.936 }, { 17.089, 14.516, 22.58 }, { 15.427, 14.336, 19.205 }, { 12.767, 11.942, 17.92 }, { 11.971, 11.553, 14.286 }, { 8.973, 10.04, 12.475 }, { 9.147, 9.339, 8.722 }, { 6.283, 8.177, 6.48 }, };

        return new_protein(num_atoms, (double *) positions);
}

struct protein *new_protein_2gb1(void)
{
        const size_t num_atoms = 56;
        const double positions[][3] = {{ -13.296, 0.028, 3.924 }, { -9.669, -0.447, 4.998 }, { -7.173, -2.314, 2.811 }, { -3.922, -3.881, 4.044 }, { -0.651, -2.752, 2.466 }, { 2.338, -5.105, 2.255 }, { 5.682, -3.321, 1.9 }, { 8.137, -5.541, 0.03 }, { 10.92, -2.963, 0.07 }, { 14.315, -4.474, -0.703 }, { 16.093, -3.026, 2.321 }, { 12.799, -2.608, 4.198 }, { 9.579, -4.606, 4.659 }, { 6.374, -3.757, 6.521 }, { 2.583, -3.604, 6.342 }, { -0.108, -1.143, 7.43 }, { -3.848, -0.651, 6.886 }, { -5.35, 2.653, 5.692 }, { -8.945, 3.892, 5.458 }, { -10.035, 4.811, 1.92 }, { -13.437, 5.248, 0.258 }, { -12.201, 3.975, -3.121 }, { -9.237, 2.226, -4.777 }, { -7.956, 5.461, -6.338 }, { -7.46, 7.449, -3.135 }, { -6.08, 4.229, -1.648 }, { -3.26, 4.204, -4.207 }, { -2.011, 7.699, -3.277 }, { -2.843, 7.366, 0.433 }, { -0.555, 4.359, 0.858 }, { 1.993, 5.797, -1.574 }, { 2.188, 8.829, 0.712 }, { 3.147, 6.462, 3.535 }, { 5.735, 4.444, 1.604 }, { 7.291, 7.781, 0.621 }, { 8.013, 8.481, 4.295 }, { 10.155, 5.328, 4.045 }, { 11.601, 5.491, 0.523 }, { 9.672, 2.625, -1.112 }, { 8.48, 3.669, -4.587 }, { 8.449, 1.498, -7.716 }, { 5.817, -1.093, -8.656 }, { 2.365, -1.045, -7.038 }, { -0.167, -3.897, -6.929 }, { -3.762, -4.276, -5.696 }, { -5.366, -7.437, -4.268 }, { -9.109, -6.81, -4.634 }, { -9.462, -10.504, -3.768 }, { -8.885, -9.951, -0.032 }, { -8.673, -6.154, 0.203 }, { -4.889, -6.162, 0.598 }, { -2.612, -3.587, -1.056 }, { 0.96, -4.492, -1.948 }, { 3.999, -2.641, -3.279 }, { 7.534, -3.776, -4.109 }, { 10.556, -1.553, -4.826 }, };

        return new_protein(num_atoms, (double *) positions);
}
