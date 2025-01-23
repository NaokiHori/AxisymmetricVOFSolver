#include <stddef.h> // size_t
#include "logger.h" // LOGGER_FAILURE
#include "domain.h" // Z_PERIODIC
#include "exchange_halo.h"

int exchange_halo_z(
    const domain_t * const domain,
    double * const * const array
) {
  if (!Z_PERIODIC) {
    LOGGER_FAILURE("z direction is not periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    array[i][     0] = array[i][nz];
    array[i][nz + 1] = array[i][ 1];
  }
  return 0;
abort:
  return 1;
}

