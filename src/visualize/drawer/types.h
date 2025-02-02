#if !defined(TYPES_H)
#define TYPES_H

#include <stdint.h>
#include <stddef.h>
#include "./vector.h"

typedef struct {
  uint8_t grayscale;
} color_t;

typedef struct {
  double depth;
  color_t color;
} pixel_t;

typedef struct {
  vector_t position;
  vector_t look_at;
} camera_t;

// triangle element which is rendered
typedef struct {
  // vertex positions
  vector_t vertices[3];
  // vertex normals (defined at each vertex)
  vector_t vertex_normals[3];
} triangle_t;

typedef struct {
  // center position
  vector_t center;
  // number of pixels
  size_t width;
  size_t height;
  // horizontal and vertical vectors of the screen horizontal and vertical axes
  vector_t local_x;
  vector_t local_y;
  // normal vector perpendicular to the screen
  vector_t normal;
} screen_t;

#endif // TYPES_H
