#if !defined(UPDATE_UX_H)
#define UPDATE_UX_H

#include "domain.h"

extern int update_ux(
    const domain_t * const domain,
    double ** const dux,
    double ** const density[2],
    double ** const ux
);

#endif // UPDATE_UX_H
