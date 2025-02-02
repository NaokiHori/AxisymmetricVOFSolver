#include <math.h>
#include "memory.h"
#include "array.h"
#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "impose_bc.h"
#include "exchange_halo.h"

static int allocate(
    const domain_t * const domain,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->ux);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->uz);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->p);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->density[0]);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->density[1]);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->txx);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->tyy);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->txz);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_field->tzz);
  return 0;
}

static int init_ux(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    for (size_t k = 0; k <= nz + 1; k++) {
      ux[i][k] = 0.;
    }
  }
  if (0 != impose_bc_ux_x(domain, ux)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (ux)");
    goto abort;
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, ux)) {
      LOGGER_FAILURE("failed to exchange halo in z (ux)");
      goto abort;
    }
  } else {
    if (0 != impose_bc_ux_z(domain, ux)) {
      LOGGER_FAILURE("failed to impose boundary condition in z (ux)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_uz(
    const domain_t * const domain,
    double ** const uz
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 0; i <= nx + 1; i++) {
    for (size_t k = 0; k <= nz + 1; k++) {
      uz[i][k] = 0.;
    }
  }
  if (0 != impose_bc_uz_x(domain, uz)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (uz)");
    goto abort;
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, uz)) {
      LOGGER_FAILURE("failed to exchange halo in z (uz)");
      goto abort;
    }
  } else {
    if (0 != impose_bc_uz_z(domain, uz)) {
      LOGGER_FAILURE("failed to impose boundary condition in z (uz)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_p(
    const domain_t * const domain,
    double ** const p
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      p[i][k] = 0.;
    }
  }
  if (0 != impose_bc_p_x(domain, p)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (p)");
    goto abort;
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, p)) {
      LOGGER_FAILURE("failed to exchange halo in z (p)");
      goto abort;
    }
  } else {
    if (0 != impose_bc_p_z(domain, p)) {
      LOGGER_FAILURE("failed to impose boundary condition in z (p)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
) {
  if (0 != allocate(domain, flow_field)) {
    LOGGER_FAILURE("failed to allocate flow field");
    goto abort;
  }
  if (0 != init_ux(domain, flow_field->ux)) {
    LOGGER_FAILURE("failed to initialise ux");
    goto abort;
  }
  if (0 != init_uz(domain, flow_field->uz)) {
    LOGGER_FAILURE("failed to initialise uz");
    goto abort;
  }
  if (0 != init_p(domain, flow_field->p)) {
    LOGGER_FAILURE("failed to initialise p");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

int flow_field_finalize(
    flow_field_t * const flow_field
) {
  destroy_array(&flow_field->ux);
  destroy_array(&flow_field->uz);
  destroy_array(&flow_field->p);
  destroy_array(&flow_field->density[0]);
  destroy_array(&flow_field->density[1]);
  destroy_array(&flow_field->txx);
  destroy_array(&flow_field->tyy);
  destroy_array(&flow_field->txz);
  destroy_array(&flow_field->tzz);
  return 0;
}

