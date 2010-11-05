#ifndef MOLECULE_H
#define MOLECULE_H
/**
 * @file molecule.h
 * @brief Interface for molecule manipulation features.
 */

#include <stddef.h>


typedef double real;

typedef real vector[3];

struct molecule {
        size_t num_atoms;  // Number of alpha carbons.
        vector *positions; // Position of each atom in $\mathbb{R}^3$.
};


extern struct molecule *new_molecule(size_t num_atoms, const vector *positions);

extern struct molecule *new_molecule_1pgb(void);

extern struct molecule *new_molecule_2gb1(void);

extern void delete_molecule(struct molecule *self);
#endif // !MOLECULE_H
