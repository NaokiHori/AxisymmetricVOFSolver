#if !defined(IMPOSE_BC_H)
#define IMPOSE_BC_H

#include "domain.h"

extern int impose_bc_ux_x(
    const domain_t * const domain,
    double * const * const ux
);

extern int impose_bc_ux_z(
    const domain_t * const domain,
    double * const * const ux
);

extern int impose_bc_uz_x(
    const domain_t * const domain,
    double * const * const uz
);

extern int impose_bc_uz_z(
    const domain_t * const domain,
    double * const * const uz
);

extern int impose_bc_p_x(
    const domain_t * const domain,
    double * const * const p
);

extern int impose_bc_p_z(
    const domain_t * const domain,
    double * const * const p
);

extern int impose_bc_vof_x(
    const domain_t * const domain,
    double * const * const vof
);

extern int impose_bc_vof_z(
    const domain_t * const domain,
    double * const * const vof
);

#endif // IMPOSE_BC_H
