#include "param.h"
#include "domain.h"
#include "exchange_halo.h"
#include "./correct_uz.h"

int correct_uz(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hzzf = domain->hzzf;
  double ** const uz = flow_field->uz;
  double ** const density = flow_field->density[1];
  double ** const psi0 = flow_solver->psi[0];
  double ** const psi1 = flow_solver->psi[1];
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[i    ][k - 1]
        + 0.5 * density[i    ][k    ];
      uz[i][k] +=
        - dt / hzzf[k] * (
            + 1. / l_density
            - 1. / DENSITY_RATIO
        ) * (
            - psi0[i    ][k - 1]
            + psi0[i    ][k    ]
        )
        - dt / hzzf[k] * (
            + 1. / DENSITY_RATIO
        ) * (
            - psi1[i    ][k - 1]
            + psi1[i    ][k    ]
        );
    }
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, uz)) {
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

