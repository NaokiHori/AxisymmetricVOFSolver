#include "logger.h"
#include "exchange_halo.h"
#include "impose_bc.h"
#include "./update_uz.h"

int update_uz(
    const domain_t * const domain,
    double ** const duz,
    double ** const density[2],
    double ** const uz
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[0][i    ][k - 1]
        + 0.5 * density[0][i    ][k    ];
      uz[i][k] = l_density * uz[i][k] + duz[i][k];
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[1][i    ][k - 1]
        + 0.5 * density[1][i    ][k    ];
      uz[i][k] /= l_density;
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

