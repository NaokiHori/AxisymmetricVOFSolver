#include "impose_bc.h"

int impose_bc_p_x(
    const domain_t * const domain,
    double * const * const p
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t k = 0; k <= nz + 1; k++) {
    p[     0][k] = p[ 1][k];
    p[nx + 1][k] = p[nx][k];
  }
  return 0;
}

int impose_bc_p_z(
    const domain_t * const domain,
    double * const * const p
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    p[i][     0] = p[i][ 1];
    p[i][nz + 1] = p[i][nz];
  }
  return 0;
}

