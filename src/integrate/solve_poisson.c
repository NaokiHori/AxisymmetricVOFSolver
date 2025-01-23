#include "param.h"
#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "rdft.h"
#include "dct.h"
#include "tdm.h"
#include "transpose.h"
#include "exchange_halo.h"
#include "./solve_poisson.h"

int solve_poisson(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  double * const buf0 = poisson_solver->buf0;
  double * const buf1 = poisson_solver->buf1;
  // assign right-hand side of Poisson equation
  {
    const double * const hxxc = domain->hxxc;
    const double * const hxxf = domain->hxxf;
    const double * const hyxc = domain->hyxc;
    const double * const hyxf = domain->hyxf;
    const double * const hzzc = domain->hzzc;
    const double * const hzzf = domain->hzzf;
    double ** const ux = flow_field->ux;
    double ** const uz = flow_field->uz;
    double ** const density = flow_field->density[1];
    double ** const psi = flow_solver->psi[0];
#pragma omp parallel for
    for (size_t i = 1; i <= nx; i++) {
      for (size_t k = 1; k <= nz; k++) {
        const double hx_xm = hxxf[i    ];
        const double hx_x0 = hxxc[i    ];
        const double hx_xp = hxxf[i + 1];
        const double hy_xm = hyxf[i    ];
        const double hy_x0 = hyxc[i    ];
        const double hy_xp = hyxf[i + 1];
        const double hz_zm = hzzf[k    ];
        const double hz_z0 = hzzc[k    ];
        const double hz_zp = hzzf[k + 1];
        const double ux_xm = ux[i    ][k    ];
        const double ux_xp = ux[i + 1][k    ];
        const double uz_zm = uz[i    ][k    ];
        const double uz_zp = uz[i    ][k + 1];
        const double den_xm = + 0.5 * density[i - 1][k    ] + 0.5 * density[i    ][k    ];
        const double den_xp = + 0.5 * density[i    ][k    ] + 0.5 * density[i + 1][k    ];
        const double den_zm = + 0.5 * density[i    ][k - 1] + 0.5 * density[i    ][k    ];
        const double den_zp = + 0.5 * density[i    ][k    ] + 0.5 * density[i    ][k + 1];
        const double dpsidx_xm = (- psi[i - 1][k    ] + psi[i    ][k    ]) / hx_xm;
        const double dpsidx_xp = (- psi[i    ][k    ] + psi[i + 1][k    ]) / hx_xp;
        const double dpsidz_zm = (- psi[i    ][k - 1] + psi[i    ][k    ]) / hz_zm;
        const double dpsidz_zp = (- psi[i    ][k    ] + psi[i    ][k + 1]) / hz_zp;
        buf0[(i - 1) * nz + (k - 1)] =
          + DENSITY_RATIO / dt / jd(hx_x0, hy_x0, hz_z0) * (
              - jdhx(hy_xm, hz_z0) * ux_xm
              + jdhx(hy_xp, hz_z0) * ux_xp
              - jdhz(hx_x0, hy_x0) * uz_zm
              + jdhz(hx_x0, hy_x0) * uz_zp
          )
          - 1. / jd(hx_x0, hy_x0, hz_z0) * (
              - jdhx(hy_xm, hz_z0) * (DENSITY_RATIO / den_xm - 1.) * dpsidx_xm
              + jdhx(hy_xp, hz_z0) * (DENSITY_RATIO / den_xp - 1.) * dpsidx_xp
              - jdhz(hx_x0, hy_x0) * (DENSITY_RATIO / den_zm - 1.) * dpsidz_zm
              + jdhz(hx_x0, hy_x0) * (DENSITY_RATIO / den_zp - 1.) * dpsidz_zp
          );
      }
    }
  }
  // project z to wave space
  if (Z_PERIODIC) {
    if (0 != rdft_exec_f(poisson_solver->rdft_plan, buf0)) {
      LOGGER_FAILURE("failed to perform RDFT");
      goto abort;
    }
  } else {
    if (0 != dct_exec_f(poisson_solver->dct_plan, buf0)) {
      LOGGER_FAILURE("failed to perform DCT2");
      goto abort;
    }
  }
  // z-align to x-align
  if (0 != transpose(nz, nx, buf0, buf1)) {
    LOGGER_FAILURE("failed to transpose array from z-aligned to x-aligned");
    goto abort;
  }
  // solve linear systems in x
  {
    tdm_plan_t * const tdm_plan = poisson_solver->tdm_plan;
    const double * const tdm_l = poisson_solver->tdm_l;
    const double * const tdm_c = poisson_solver->tdm_c;
    const double * const tdm_u = poisson_solver->tdm_u;
    const double * const wavenumbers = poisson_solver->wavenumbers;
    if (0 != tdm_solve(tdm_plan, tdm_l, tdm_c, tdm_u, wavenumbers, buf1)) {
      LOGGER_FAILURE("failed to solve tri-diagonal matrix");
      goto abort;
    }
  }
  // x-align to z-align
  if (0 != transpose(nx, nz, buf1, buf0)) {
    LOGGER_FAILURE("failed to transpose array from x-aligned to z-aligned");
    goto abort;
  }
  // project z to physical space
  if (Z_PERIODIC) {
    if (0 != rdft_exec_b(poisson_solver->rdft_plan, buf0)) {
      LOGGER_FAILURE("failed to perform IRDFT");
      goto abort;
    }
  } else {
    if (0 != dct_exec_b(poisson_solver->dct_plan, buf0)) {
      LOGGER_FAILURE("failed to perform DCT3");
      goto abort;
    }
  }
  {
    double ** const psi = flow_solver->psi[1];
#pragma omp parallel for
    for (size_t i = 1; i <= nx; i++) {
      for (size_t k = 1; k <= nz; k++) {
        psi[i][k] = buf0[(i - 1) * nz + (k - 1)];
      }
    }
    if (Z_PERIODIC) {
      if (0 != exchange_halo_z(domain, psi)) {
        LOGGER_FAILURE("failed to exchange halo in z");
        goto abort;
      }
    } else {
      for (size_t i = 1; i <= nx; i++) {
        psi[i][     0] = psi[i][ 1];
        psi[i][nz + 1] = psi[i][nz];
      }
    }
    for (size_t k = 0; k <= nz + 1; k++) {
      psi[     0][k] = psi[ 1][k];
      psi[nx + 1][k] = psi[nx][k];
    }
  }
  return 0;
abort:
  return 1;
}

