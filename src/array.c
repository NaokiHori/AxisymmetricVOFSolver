#include "memory.h"
#include "array.h"

int prepare_array(
    const size_t nitems[NDIMS],
    double *** const array
) {
  double * const buffer = memory_alloc(nitems[0] * nitems[1], sizeof(double));
  *array = memory_alloc(nitems[0], sizeof(double *));
  for (size_t i = 0; i < nitems[0]; i++) {
    (*array)[i] = buffer + nitems[1] * i;
  }
  return 0;
}

int destroy_array(
    double *** const array
) {
  memory_free((*array)[0]);
  memory_free(*array);
  *array = NULL;
  return 0;
}

