#if !defined(RENDER_H)
#define RENDER_H

#include "./types.h"

extern int render_triangle(
    const camera_t * const camera,
    const vector_t * const light,
    const screen_t * const screen,
    const triangle_t * const triangle,
    pixel_t * const canvas
);

#endif // RENDER_H
