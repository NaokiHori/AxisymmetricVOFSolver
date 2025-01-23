#include "impose_bc.h"

int impose_bc_vof_x(
    const domain_t * const domain,
    double * const * const vof
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t k = 0; k <= nz + 1; k++) {
    vof[     0][k] = vof[ 1][k];
    vof[nx + 1][k] = vof[nx][k];
  }
  return 0;
}

int impose_bc_vof_z(
    const domain_t * const domain,
    double * const * const vof
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    vof[i][     0] = 0.;
    vof[i][nz + 1] = 1.;
  }
  return 0;
}

