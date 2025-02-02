#include <math.h>
#include "array.h"
#include "impose_bc.h"
#include "exchange_halo.h"
#include "interface_field.h"

int interface_field_init(
    const domain_t * const domain,
    interface_field_t * const interface_field
) {
  const double lz = domain->lz;
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const xc = domain->xc;
  const double * const zc = domain->zc;
  const double * const hxxc = domain->hxxc;
  const double * const hzzc = domain->hzzc;
  double *** const vof = &interface_field->vof;
  double * const vof_beta = &interface_field->vof_beta;
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, vof);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->vof_ux);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->vof_uz);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->corner_nx);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->corner_nz);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->center_nx);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->center_nz);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->curv);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->force_x);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &interface_field->force_z);
  *vof_beta = 1.;
  const double bath_height = 0.25 * lz;
  const double cz = 0.75 * lz;
  const double radius = 0.5;
  // water: 0
  // air: 1
  // filled with water if
  //   z < bath_height
  //   distance from (0, cz) is less than radius
  for (size_t i = 1; i <= nx; i++) {
    const double x = xc[i];
    for (size_t k = 1; k <= nz; k++) {
      const double z = zc[k];
      const double bath_distance = z - bath_height;
      const double drop_distance = sqrt(pow(x, 2.) + pow(z - cz, 2.)) - radius;
      const double distance = fmin(bath_distance, drop_distance) / fmin(hxxc[i], hzzc[k]);
      (*vof)[i][k] = 0.5 * (1. + tanh(*vof_beta * distance));
    }
  }
  impose_bc_vof_x(domain, *vof);
  if (Z_PERIODIC) {
    exchange_halo_z(domain, *vof);
  } else {
    impose_bc_vof_z(domain, *vof);
  }
  return 0;
}

int interface_field_finalize(
    interface_field_t * const interface_field
) {
  destroy_array(&interface_field->vof);
  destroy_array(&interface_field->vof_ux);
  destroy_array(&interface_field->vof_uz);
  destroy_array(&interface_field->corner_nx);
  destroy_array(&interface_field->corner_nz);
  destroy_array(&interface_field->center_nx);
  destroy_array(&interface_field->center_nz);
  destroy_array(&interface_field->curv);
  destroy_array(&interface_field->force_x);
  destroy_array(&interface_field->force_z);
  return 0;
}

