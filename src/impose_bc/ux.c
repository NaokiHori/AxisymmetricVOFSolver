#include <stddef.h> // size_t
#include "logger.h"
#include "domain.h"
#include "impose_bc.h"

int impose_bc_ux_x(
    const domain_t * const domain,
    double * const * const ux
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t k = 0; k <= nz + 1; k++) {
    ux[     0][k] = 0.; // unused
    ux[     1][k] = 0.;
    ux[nx + 1][k] = 0.;
  }
  return 0;
}

int impose_bc_ux_z(
    const domain_t * const domain,
    double * const * const ux
) {
  if (Z_PERIODIC) {
    LOGGER_FAILURE("z direction is periodic");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    ux[i][     0] = 0.;
    ux[i][nz + 1] = 0.;
  }
  return 0;
abort:
  return 1;
}

