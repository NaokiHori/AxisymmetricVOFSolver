#if !defined(UPDATE_UZ_H)
#define UPDATE_UZ_H

#include "domain.h"

extern int update_uz(
    const domain_t * const domain,
    double ** const duz,
    double ** const density[2],
    double ** const uz
);

#endif // UPDATE_UZ_H
