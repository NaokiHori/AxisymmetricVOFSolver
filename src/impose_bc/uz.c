#include <stddef.h> // size_t
#include "logger.h"
#include "domain.h"
#include "impose_bc.h"

int impose_bc_uz_x(
    const domain_t * const domain,
    double * const * const uz
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t k = 0; k <= nz + 1; k++) {
    uz[     0][k] = uz[ 1][k];
    uz[nx + 1][k] = 0.;
  }
  return 0;
}

int impose_bc_uz_z(
    const domain_t * const domain,
    double * const * const uz
) {
  if (Z_PERIODIC) {
    LOGGER_FAILURE("z direction is periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    uz[i][     0] = 0.; // unused
    uz[i][     1] = 0.;
    uz[i][nz + 1] = 0.;
  }
  return 0;
abort:
  return 1;
}

