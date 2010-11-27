#ifndef POTENTIAL_H
#define POTENTIAL_H

struct protein;

struct contact_map;


extern double potential(const struct protein *p,
                        const struct contact_map *native_map,
                        double a);

#endif // POTENTIAL_H
