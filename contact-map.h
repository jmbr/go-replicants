#ifndef CONTACT_MAP_H
#define CONTACT_MAP_H   1
/**
 * @file contact_map.h
 */

#include <stdio.h>
#include <stddef.h>


struct protein;

struct contact_map;


extern struct contact_map *new_contact_map(const struct protein *protein,
                                           double d_max);

extern void delete_contact_map(struct contact_map *self);

extern size_t contact_map_get_num_contacts(const struct contact_map *self);

extern double contact_map_get_d_max(const struct contact_map *self);

extern double contact_map_get_distance(const struct contact_map *self,
                                       size_t i, size_t j);

extern int contact_map_plot(const struct contact_map *self, FILE *gnuplot);
#endif // !CONTACT_MAP_H
