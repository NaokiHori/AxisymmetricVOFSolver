#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include "dct.h" // dct_plan_t
#include "rdft.h" // rdft_plan_t
#include "tdm.h" // tdm_plan_t
#include "domain.h" // domain_t

// variables used to solve Poisson equations
typedef struct {
  // buffers to store intermediate data
  double * buf0;
  double * buf1;
  // x direction: tdm-related things
  tdm_plan_t * tdm_plan;
  // z direction: dft-related things
  rdft_plan_t * rdft_plan;
  dct_plan_t * dct_plan;
  double * wavenumbers;
  double * tdm_l;
  double * tdm_c;
  double * tdm_u;
} poisson_solver_t;

typedef struct {
  double ** psi[2];
  double ** dux;
  double ** duz;
  poisson_solver_t poisson_solver;
} flow_solver_t;

extern int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
);

extern int flow_solver_finalize(
    flow_solver_t * const flow_solver
);

#endif // FLOW_SOLVER_H
