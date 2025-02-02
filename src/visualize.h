#if !defined(VISUALIZE_H)
#define VISUALIZE_H

#include <stddef.h>
#include <stdbool.h>

typedef struct {
  double x;
  double y;
  double nx;
  double ny;
} point_t;

typedef struct {
  bool is_closed;
  size_t npoints;
  point_t * points;
} contour_t;

extern int visualize(
    const char file_name[],
    const size_t nx,
    const size_t nz,
    const double * const x,
    const double * const z,
    const double * const * const vof
);

#endif // VISUALIZE_H
