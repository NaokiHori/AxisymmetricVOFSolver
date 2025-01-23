#if !defined(CORRECT_UX_H)
#define CORRECT_UX_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"

extern int correct_ux(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
);

#endif // CORRECT_UX_H
