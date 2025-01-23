#if !defined(SAVE_H)
#define SAVE_H

#include <stddef.h> // size_t
#include "domain.h" // domain_t
#include "flow_field.h" // flow_field_t
#include "interface_field.h" // interface_field_t

extern int save(
    const size_t id,
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field
);

#endif // SAVE_H
