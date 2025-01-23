#if !defined(COMPUTE_DUZ_H)
#define COMPUTE_DUZ_H

#include "domain.h"
#include "flow_field.h"
#include "interface_field.h"

int compute_duz(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field,
    const double dt,
    double ** const force_z,
    double ** const duz
);

#endif // COMPUTE_DUZ_H
