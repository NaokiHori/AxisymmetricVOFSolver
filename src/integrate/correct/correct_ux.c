#include "param.h"
#include "domain.h"
#include "exchange_halo.h"
#include "./correct_ux.h"

int correct_ux(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxf = domain->hxxf;
  double ** const ux = flow_field->ux;
  double ** const density = flow_field->density[1];
  double ** const psi0 = flow_solver->psi[0];
  double ** const psi1 = flow_solver->psi[1];
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[i - 1][k    ]
        + 0.5 * density[i    ][k    ];
      ux[i][k] +=
        - dt / hxxf[i] * (
            + 1. / l_density
            - 1. / DENSITY_RATIO
        ) * (
            - psi0[i - 1][k    ]
            + psi0[i    ][k    ]
        )
        - dt / hxxf[i] * (
            + 1. / DENSITY_RATIO
        ) * (
            - psi1[i - 1][k    ]
            + psi1[i    ][k    ]
        );
    }
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, ux)) {
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

