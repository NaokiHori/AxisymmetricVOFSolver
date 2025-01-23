#if !defined(COMPUTE_PHYSICAL_PROPERTY_H)
#define COMPUTE_PHYSICAL_PROPERTY_H

#include "domain.h"
#include "interface_field.h"

extern int compute_density(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    double ** const density
);

#endif // COMPUTE_PHYSICAL_PROPERTY_H
