#if !defined(INTERFACE_FIELD_H)
#define INTERFACE_FIELD_H

#include "domain.h" // domain_t

typedef struct {
  double ** vof;
  double ** vof_ux;
  double ** vof_uz;
  double ** corner_nx;
  double ** corner_nz;
  double ** center_nx;
  double ** center_nz;
  double ** curv;
  double ** force_x;
  double ** force_z;
  double vof_beta;
} interface_field_t;

extern int interface_field_init(
    const domain_t * const domain,
    interface_field_t * const interface_field
);

extern int interface_field_finalize(
    interface_field_t * const interface_field
);

#endif // INTERFACE_FIELD_H
