#include <stddef.h>
#include "param.h"
#include "./compute_density.h"

static double clamp(
    double vof
) {
  vof = vof < 0. ? 0. : vof;
  vof = 1. < vof ? 1. : vof;
  return vof;
}

int compute_density(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    double ** const density
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  double ** const vof = interface_field->vof;
#pragma omp parallel for
  for (size_t i = 0; i <= nx + 1; i++) {
    for (size_t k = 0; k <= nz + 1; k++) {
      const double clamped_vof = clamp(vof[i][k]);
      density[i][k] = 1. + (DENSITY_RATIO - 1.) * clamped_vof;
    }
  }
  return 0;
}

