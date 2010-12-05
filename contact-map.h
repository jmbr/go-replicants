#ifndef CONTACT_MAP_H
#define CONTACT_MAP_H

struct protein;

struct contact_map;


extern struct contact_map *new_contact_map(const struct protein *protein,
                                           double d_max);
extern void delete_contact_map(struct contact_map *self);

extern size_t contact_map_get_num_contacts(const struct contact_map *self);
extern double contact_map_get_d_max(const struct contact_map *self);

extern double contact_map_get_distance(const struct contact_map *self,
                                       size_t i, size_t j);

extern void contact_map_plot(const struct contact_map *self, FILE *gnuplot,
                             const char *title_format, ...)
        __attribute__ ((format (printf, 3, 4)));

extern int contact_map_diff(const struct contact_map *c1, const struct contact_map *c2);

#endif // !CONTACT_MAP_H
