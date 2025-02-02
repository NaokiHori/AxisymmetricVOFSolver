#if !defined(DRAWER_H)
#define DRAWER_H

#include <stddef.h>
#include "../visualize.h"

extern int draw(
    const char file_name[],
    const size_t ncontours,
    const contour_t * const contours
);

#endif // DRAWER_H
