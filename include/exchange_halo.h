#if !defined(EXCHANGE_HALO_H)
#define EXCHANGE_HALO_H

#include "domain.h"

extern int exchange_halo_z(
    const domain_t * const domain,
    double * const * const array
);

#endif // EXCHANGE_HALO_H
