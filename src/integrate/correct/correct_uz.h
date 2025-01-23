#if !defined(CORRECT_UZ_H)
#define CORRECT_UZ_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"

extern int correct_uz(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

#endif // CORRECT_UZ_H
