#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "./correct.h"
#include "./correct/correct_ux.h"
#include "./correct/correct_uz.h"

int correct(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  if (0 != correct_ux(domain, flow_field, flow_solver, dt)) {
    LOGGER_FAILURE("failed to correct uz");
    goto abort;
  }
  if (0 != correct_uz(domain, flow_field, flow_solver, dt)) {
    LOGGER_FAILURE("failed to correct uz");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to correct velocity field");
  return 1;
}

