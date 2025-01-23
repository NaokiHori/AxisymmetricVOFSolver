#if !defined(COMPUTE_SHEAR_STRESS_TENSOR_H)
#define COMPUTE_SHEAR_STRESS_TENSOR_H

#include "domain.h"
#include "flow_field.h"
#include "interface_field.h"

extern int compute_shear_stress_tensor(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    flow_field_t * const flow_field
);

#endif // COMPUTE_SHEAR_STRESS_TENSOR_H
