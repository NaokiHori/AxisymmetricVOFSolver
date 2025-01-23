#include "logger.h"
#include "exchange_halo.h"
#include "impose_bc.h"
#include "./update_ux.h"

int update_ux(
    const domain_t * const domain,
    double ** const dux,
    double ** const density[2],
    double ** const ux
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[0][i - 1][k    ]
        + 0.5 * density[0][i    ][k    ];
      ux[i][k] = l_density * ux[i][k] + dux[i][k];
    }
  }
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[1][i - 1][k    ]
        + 0.5 * density[1][i    ][k    ];
      ux[i][k] /= l_density;
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

