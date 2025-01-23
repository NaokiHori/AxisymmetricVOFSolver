#include <math.h>
#include "memory.h"
#include "domain.h"

const size_t ux_imin = 2;
const size_t uz_kmin = Z_PERIODIC ? 1 : 2;

int domain_init(
    domain_t * const domain
) {
  double * const lx = &domain->lx;
  double * const lz = &domain->lz;
  size_t * const nx = &domain->nx;
  size_t * const nz = &domain->nz;
  double ** const xf = &domain->xf;
  double ** const xc = &domain->xc;
  double ** const zf = &domain->zf;
  double ** const zc = &domain->zc;
  double ** const hxxf = &domain->hxxf;
  double ** const hxxc = &domain->hxxc;
  double ** const hyxf = &domain->hyxf;
  double ** const hyxc = &domain->hyxc;
  double ** const hzzf = &domain->hzzf;
  double ** const hzzc = &domain->hzzc;
  *lx = 2.;
  *lz = 8.;
  *nx = 32;
  *nz = 128;
  *xf = memory_alloc(*nx + 2, sizeof(double));
  *xc = memory_alloc(*nx + 2, sizeof(double));
  *zf = memory_alloc(*nz + 2, sizeof(double));
  *zc = memory_alloc(*nz + 2, sizeof(double));
  *hxxf = memory_alloc(*nx + 2, sizeof(double));
  *hxxc = memory_alloc(*nx + 2, sizeof(double));
  *hyxf = memory_alloc(*nx + 2, sizeof(double));
  *hyxc = memory_alloc(*nx + 2, sizeof(double));
  *hzzf = memory_alloc(*nz + 2, sizeof(double));
  *hzzc = memory_alloc(*nz + 2, sizeof(double));
  // positions
  for (size_t i = 0; i <= *nx + 1; i++) {
    if (0 == i) {
      (*xf)[i] = nan("");
    } else {
      (*xf)[i] = 1. * (i - 1) * (*lx / *nx);
    }
  }
  for (size_t i = 0; i <= *nx + 1; i++) {
    if (0 == i) {
      (*xc)[i] = 0.;
    } else if (*nx + 1 == i) {
      (*xc)[i] = *lx;
    } else {
      (*xc)[i] =
        + 0.5 * (*xf)[i    ]
        + 0.5 * (*xf)[i + 1];
    }
  }
  for (size_t k = 0; k <= *nz + 1; k++) {
    if (0 == k) {
      (*zf)[k] = nan("");
    } else {
      (*zf)[k] = 1. * (k - 1) * (*lz / *nz);
    }
  }
  for (size_t k = 0; k <= *nz + 1; k++) {
    if (0 == k) {
      (*zc)[k] = 0.;
    } else if (*nz + 1 == k) {
      (*zc)[k] = *lz;
    } else {
      (*zc)[k] =
        + 0.5 * (*zf)[k    ]
        + 0.5 * (*zf)[k + 1];
    }
  }
  // scale factors
  for (size_t i = 0; i <= *nx + 1; i++) {
    if (0 == i) {
      (*hxxf)[i] = nan("");
    } else {
      (*hxxf)[i] = - (*xc)[i - 1]
                   + (*xc)[i    ];
    }
  }
  for (size_t i = 0; i <= *nx + 1; i++) {
    if (0 == i) {
      (*hxxc)[i] = - (*xc)[0]
                   + (*xc)[1];
    } else if (*nx + 1 == i) {
      (*hxxc)[i] = - (*xc)[*nx    ]
                   + (*xc)[*nx + 1];
    } else {
      (*hxxc)[i] = - (*xf)[i    ]
                   + (*xf)[i + 1];
    }
  }
  for (size_t i = 0; i <= *nx + 1; i++) {
    if (0 == i) {
      (*hyxf)[i] = nan("");
    } else {
      (*hyxf)[i] = (*xf)[i];
    }
  }
  for (size_t i = 0; i <= *nx + 1; i++) {
    (*hyxc)[i] = (*xc)[i];
  }
  if (Z_PERIODIC) {
    for (size_t k = 0; k <= *nz + 1; k++) {
      (*hzzf)[k] = *lz / *nz;
    }
    for (size_t k = 0; k <= *nz + 1; k++) {
      (*hzzc)[k] = *lz / *nz;
    }
  } else {
    for (size_t k = 0; k <= *nz + 1; k++) {
      if (0 == k) {
        (*hzzf)[k] = nan("");
      } else if (1 == k || *nz + 1 == k) {
        (*hzzf)[k] = 0.5 * *lz / *nz;
      } else {
        (*hzzf)[k] = *lz / *nz;
      }
    }
    for (size_t k = 0; k <= *nz + 1; k++) {
      if (0 == k || *nz + 1 == k) {
        (*hzzc)[k] = 0.5 * *lz / *nz;
      } else {
        (*hzzc)[k] = *lz / *nz;
      }
    }
  }
  return 0;
}

int domain_finalize(
    domain_t * const domain
) {
  memory_free(domain->xf);
  memory_free(domain->xc);
  memory_free(domain->zf);
  memory_free(domain->zc);
  memory_free(domain->hxxf);
  memory_free(domain->hxxc);
  memory_free(domain->hyxf);
  memory_free(domain->hyxc);
  memory_free(domain->hzzf);
  memory_free(domain->hzzc);
  return 0;
}

double jd(
    const double hx,
    const double hy,
    const double hz
) {
  return hx * hy * hz;
}

double jdhx(
    const double hy,
    const double hz
) {
  return hy * hz;
}

double jdhz(
    const double hx,
    const double hy
) {
  return hx * hy;
}

