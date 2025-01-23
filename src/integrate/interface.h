#if !defined(INTEGRATE_INTERFACE_H)
#define INTEGRATE_INTERFACE_H

#include "flow_field.h"
#include "interface_field.h"

extern int update_interface(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    interface_field_t * const interface_field,
    const double dt
);

#endif // INTEGRATE_INTERFACE_H
