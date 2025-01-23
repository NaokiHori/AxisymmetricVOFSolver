#include "logger.h"
#include "interface_field.h"
#include "./predict.h"
#include "./predict/compute_shear_stress_tensor.h"
#include "./predict/compute_dux.h"
#include "./predict/compute_duz.h"
#include "./predict/update_ux.h"
#include "./predict/update_uz.h"

int predict(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    interface_field_t * const interface_field,
    const double dt
) {
  double ** const force_x = interface_field->force_x;
  double ** const force_z = interface_field->force_z;
  double ** const dux = flow_solver->dux;
  double ** const duz = flow_solver->duz;
  if (0 != compute_shear_stress_tensor(domain, interface_field, flow_field)) {
    LOGGER_FAILURE("failed to compute shear-stress tensor");
    goto abort;
  }
  if (0 != compute_dux(domain, flow_field, interface_field, dt, force_x, dux)) {
    LOGGER_FAILURE("failed to find dux");
    goto abort;
  }
  if (0 != compute_duz(domain, flow_field, interface_field, dt, force_z, duz)) {
    LOGGER_FAILURE("failed to find duz");
    goto abort;
  }
  if (0 != update_ux(domain, dux, flow_field->density, flow_field->ux)) {
    LOGGER_FAILURE("failed to update ux");
    goto abort;
  }
  if (0 != update_uz(domain, duz, flow_field->density, flow_field->uz)) {
    LOGGER_FAILURE("failed to update uz");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to predict flow field");
  return 1;
}

