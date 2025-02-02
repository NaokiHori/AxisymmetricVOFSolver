#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "./vector.h"

vector_t normalize(
    const vector_t v
) {
  const double norm = sqrt(
      + pow(v.x, 2)
      + pow(v.y, 2)
      + pow(v.z, 2)
  );
  if (norm < DBL_EPSILON) {
    puts("trying to normalize a zero-norm vector");
    exit(EXIT_FAILURE);
  }
  return (vector_t){
    .x = v.x / norm,
    .y = v.y / norm,
    .z = v.z / norm,
  };
}

vector_t add(
    const vector_t v0,
    const vector_t v1
) {
  return (vector_t){
    .x = v0.x + v1.x,
    .y = v0.y + v1.y,
    .z = v0.z + v1.z,
  };
}

vector_t sub(
    const vector_t v0,
    const vector_t v1
) {
  return (vector_t){
    .x = v0.x - v1.x,
    .y = v0.y - v1.y,
    .z = v0.z - v1.z,
  };
}

vector_t mul(
    const double c,
    const vector_t v
) {
  return (vector_t){
    .x = c * v.x,
    .y = c * v.y,
    .z = c * v.z,
  };
}

double inner_product(
    const vector_t v0,
    const vector_t v1
) {
  return
    + v0.x * v1.x
    + v0.y * v1.y
    + v0.z * v1.z;
}

vector_t outer_product(
    const vector_t v0,
    const vector_t v1
) {
  return (vector_t){
    .x = v0.y * v1.z - v0.z * v1.y,
    .y = v0.z * v1.x - v0.x * v1.z,
    .z = v0.x * v1.y - v0.y * v1.x,
  };
}

