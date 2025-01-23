#include <stdio.h>
#include <math.h> // fmin, fmax, fabs, pow
#include "logger.h"
#include "param.h"
#include "domain.h"
#include "flow_field.h"
#include "./decide_dt.h"

static int decide_dt_adv(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxf = domain->hxxf;
  const double * const hzzf = domain->hzzf;
  const double small = 1.e-8;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  *dt = 1.;
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      *dt = fmin(*dt, hxxf[i] / fmax(small, fabs(ux[i][k])));
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      *dt = fmin(*dt, hzzf[k] / fmax(small, fabs(uz[i][k])));
    }
  }
  // safety factor
  *dt *= 0.125;
  return 0;
}

static int decide_dt_dif(
    const domain_t * const domain,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  double min_hx = domain->lx;
  double min_hz = domain->lz;
  const double * const hxxc = domain->hxxc;
  const double * const hzzc = domain->hzzc;
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      min_hx = fmin(min_hx, hxxc[i]);
    }
  }
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      min_hz = fmin(min_hz, hzzc[k]);
    }
  }
  *dt = (0.5 / NDIMS) * (Re * fmin(1., DENSITY_RATIO / VISCOSITY_RATIO)) * pow(fmin(min_hx, min_hz), 2.);
  // safety factor
  *dt *= 0.5;
  return 0;
}

static int decide_dt_int(
    const domain_t * const domain,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  double min_hx = domain->lx;
  double min_hz = domain->lz;
  const double * const hxxc = domain->hxxc;
  const double * const hzzc = domain->hzzc;
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      min_hx = fmin(min_hx, hxxc[i]);
    }
  }
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      min_hz = fmin(min_hz, hzzc[k]);
    }
  }
  const double pi = 3.14159265358979323846;
  *dt = sqrt(We * 0.25 / pi * pow(fmin(min_hx, min_hz), 3.));
  // safety factor
  *dt *= 0.5;
  return 0;
}

int decide_dt(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  double dt_adv = 0.;
  double dt_dif = 0.;
  double dt_int = 0.;
  if (0 != decide_dt_adv(domain, flow_field, &dt_adv)) {
    LOGGER_FAILURE("failed to find advective time-step constraint");
    goto abort;
  }
  if (0 != decide_dt_dif(domain, &dt_dif)) {
    LOGGER_FAILURE("failed to find diffusive time-step constraint");
    goto abort;
  }
  if (0 != decide_dt_int(domain, &dt_int)) {
    LOGGER_FAILURE("failed to find intrefacial time-step constraint");
    goto abort;
  }
  *dt = fmin(dt_adv, fmin(dt_dif, dt_int));
  return 0;
abort:
  LOGGER_FAILURE("failed to find time-step size");
  return 1;
}

