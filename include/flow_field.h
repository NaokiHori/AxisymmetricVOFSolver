#if !defined(FLOW_FIELD_H)
#define FLOW_FIELD_H

#include "domain.h"

typedef struct {
  double ** ux;
  double ** uz;
  double ** p;
  double ** density[2];
  double ** txx;
  double ** tyy;
  double ** txz;
  double ** tzz;
} flow_field_t;

extern int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
);

extern int flow_field_finalize(
    flow_field_t * const flow_field
);

#endif // FLOW_FIELD_H
