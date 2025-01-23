#include "domain.h"
#include "logger.h"
#include "integrate.h"
#include "./integrate/compute_density.h"
#include "./integrate/decide_dt.h"
#include "./integrate/interface.h"
#include "./integrate/predict.h"
#include "./integrate/solve_poisson.h"
#include "./integrate/correct.h"
#include "./integrate/update_pressure.h"

int integrate(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    interface_field_t * const interface_field,
    double * const dt
) {
  if (0 != compute_density(domain, interface_field, flow_field->density[0])) {
    LOGGER_FAILURE("failed to compute density");
    goto abort;
  }
  if (0 != decide_dt(domain, flow_field, dt)) {
    LOGGER_FAILURE("failed to find time-step size");
    goto abort;
  }
  if (0 != update_interface(domain, flow_field, interface_field, *dt)) {
    LOGGER_FAILURE("failed to update interface field");
    goto abort;
  }
  if (0 != compute_density(domain, interface_field, flow_field->density[1])) {
    LOGGER_FAILURE("failed to compute density");
    goto abort;
  }
  if (0 != predict(domain, flow_field, flow_solver, interface_field, *dt)) {
    LOGGER_FAILURE("failed to predict flow field");
    goto abort;
  }
  if (0 != solve_poisson(domain, flow_field, flow_solver, *dt)) {
    LOGGER_FAILURE("failed to solve Poisson equation to find scalar potential");
    goto abort;
  }
  if (0 != correct(domain, flow_field, flow_solver, *dt)) {
    LOGGER_FAILURE("failed to enforce incompressibility");
    goto abort;
  }
  if (0 != update_pressure(domain, flow_field, flow_solver)) {
    LOGGER_FAILURE("failed to update pressure field");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to update flow field");
  return 1;
}

