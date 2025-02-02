#if !defined(OUTPUT_H)
#define OUTPUT_H

#include <stddef.h>
#include "./types.h"

extern int output(
    const char file_name[],
    const size_t width,
    const size_t height,
    pixel_t * const canvas
);

#endif // OUTPUT_H
