#include <math.h>
#include "memory.h"
#include "logger.h"
#include "array.h"
#include "domain.h"
#include "flow_solver.h"
#include "rdft.h"
#include "dct.h"
#include "tdm.h"

static int init_x_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hxxf = domain->hxxf;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  tdm_plan_t ** const tdm_plan = &poisson_solver->tdm_plan;
  if (0 != tdm_init_plan(nx, nz, false, tdm_plan)) {
    LOGGER_FAILURE("failed to initialize TDM solver");
    goto abort;
  }
  double ** const tdm_l = &poisson_solver->tdm_l;
  double ** const tdm_c = &poisson_solver->tdm_c;
  double ** const tdm_u = &poisson_solver->tdm_u;
  *tdm_l = memory_alloc(nx, sizeof(double));
  *tdm_c = memory_alloc(nx, sizeof(double));
  *tdm_u = memory_alloc(nx, sizeof(double));
  for (size_t i = 1; i <= nx; i++) {
    double * const l = (*tdm_l) + i - 1;
    double * const c = (*tdm_c) + i - 1;
    double * const u = (*tdm_u) + i - 1;
    // NOTE: assuming z direction being homogeneous
    //       taking k = 1 as a representative position
    const size_t k = 1;
    const double hx_xm = hxxf[i    ];
    const double hx_x0 = hxxc[i    ];
    const double hx_xp = hxxf[i + 1];
    const double hy_xm = hyxf[i    ];
    const double hy_x0 = hyxc[i    ];
    const double hy_xp = hyxf[i + 1];
    const double hz = hzzc[k];
    *l = 1. / jd(hx_x0, hy_x0, hz) * jdhx(hy_xm, hz) / hx_xm;
    *u = 1. / jd(hx_x0, hy_x0, hz) * jdhx(hy_xp, hz) / hx_xp;
    // impose Neumann boundary condition dp/dx = 0
    if (1 == i) {
      *c = - 1. * *u;
    } else if (nx == i) {
      *c = - 1. * *l;
    } else {
      *c = - 1. * *l - 1. * *u;
    }
  }
  return 0;
abort:
  return 1;
}

static int init_z_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const double lz = domain->lz;
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double dz = lz / nz;
  double ** const wavenumbers = &poisson_solver->wavenumbers;
  *wavenumbers = memory_alloc(nz, sizeof(double));
  if (Z_PERIODIC) {
    rdft_plan_t ** const rdft_plan = &poisson_solver->rdft_plan;
    if (0 != rdft_init_plan(nz, nx, rdft_plan)) {
      LOGGER_FAILURE("failed to initialize RDFT solver");
      goto abort;
    }
  } else {
    dct_plan_t ** const dct_plan = &poisson_solver->dct_plan;
    if (0 != dct_init_plan(nz, nx, dct_plan)) {
      LOGGER_FAILURE("failed to initialize DCT solver");
      goto abort;
    }
  }
  const double dft_norm = Z_PERIODIC ? 1. * nz : 2. * nz;
  for (size_t k = 0; k < nz; k++) {
    const double pi = 3.1415926535897932385;
    (*wavenumbers)[k] = - pow(2. / dz * sin(pi * k / dft_norm), 2.);
  }
  return 0;
abort:
  return 1;
}

int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  // auxiliary buffers
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_solver->psi[0]);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_solver->psi[1]);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_solver->dux);
  prepare_array((size_t [NDIMS]){nx + 2, nz + 2}, &flow_solver->duz);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  double ** const buf0 = &poisson_solver->buf0;
  double ** const buf1 = &poisson_solver->buf1;
  *buf0 = memory_alloc(nx * nz, sizeof(double));
  *buf1 = memory_alloc(nx * nz, sizeof(double));
  // x direction: tdm-related things
  if (0 != init_x_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialize dft part of poisson solver");
    goto abort;
  }
  // z direction: dft-related things
  if (0 != init_z_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialize tdm part of poisson solver");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to initialize flow solver");
  return 1;
}

int flow_solver_finalize(
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  destroy_array(&flow_solver->psi[0]);
  destroy_array(&flow_solver->psi[1]);
  destroy_array(&flow_solver->dux);
  destroy_array(&flow_solver->duz);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  memory_free(poisson_solver->buf0);
  memory_free(poisson_solver->buf1);
  if (Z_PERIODIC) {
    rdft_destroy_plan(&poisson_solver->rdft_plan);
  } else {
    dct_destroy_plan(&poisson_solver->dct_plan);
  }
  tdm_destroy_plan(&poisson_solver->tdm_plan);
  memory_free(poisson_solver->wavenumbers);
  memory_free(poisson_solver->tdm_l);
  memory_free(poisson_solver->tdm_c);
  memory_free(poisson_solver->tdm_u);
  return 0;
}

