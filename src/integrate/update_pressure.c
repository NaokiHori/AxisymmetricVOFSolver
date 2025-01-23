#include "logger.h"
#include "exchange_halo.h"
#include "./update_pressure.h"

int update_pressure(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  {
    double ** const psi = flow_solver->psi[1];
    double ** const p = flow_field->p;
#pragma omp parallel for
    for (size_t i = 1; i <= nx; i++) {
      for (size_t k = 1; k <= nz; k++) {
        p[i][k] += psi[i][k];
      }
    }
    if (Z_PERIODIC) {
      if (0 != exchange_halo_z(domain, p)) {
        LOGGER_FAILURE("failed to exchange halo in z");
        goto abort;
      }
    }
  }
  {
    double ** const psi0 = flow_solver->psi[0];
    double ** const psi1 = flow_solver->psi[1];
#pragma omp parallel for
    for (size_t i = 0; i <= nx + 1; i++) {
      for (size_t k = 0; k <= nz + 1; k++) {
        psi0[i][k] = psi1[i][k];
      }
    }
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to update pressure field");
  return 1;
}

