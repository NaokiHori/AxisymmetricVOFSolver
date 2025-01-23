#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h> // size_t
#include <stdbool.h> // true, false

#define Z_PERIODIC false

extern const size_t ux_imin;
extern const size_t uz_kmin;

typedef struct {
  double lx;
  double lz;
  size_t nx;
  size_t nz;
  double * xf;
  double * xc;
  double * zf;
  double * zc;
  double * hxxf;
  double * hxxc;
  double * hyxf;
  double * hyxc;
  double * hzzf;
  double * hzzc;
} domain_t;

extern int domain_init(
    domain_t * const domain
);

extern int domain_finalize(
    domain_t * const domain
);

extern double jd(
    const double hx,
    const double hy,
    const double hz
);

extern double jdhx(
    const double hy,
    const double hz
);

extern double jdhz(
    const double hx,
    const double hy
);

#endif // DOMAIN_H
