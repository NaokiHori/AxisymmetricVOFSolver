#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include "logger.h"
#include "domain.h"
#include "monitor.h"

#define ROOT_DIRECTORY "output/log/"

static int output(
    const size_t step,
    const double time,
    const char file_name[],
    const size_t nitems,
    const double * const quantities
) {
  errno = 0;
  FILE * const fp = fopen(file_name, "a");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  fprintf(fp, "%10zu % .15e ", step, time);
  for (size_t n = 0; n < nitems; n++) {
    fprintf(fp, "% .15e%c", quantities[n], nitems - 1 == n ? '\n' : ' ');
  }
  fclose(fp);
  return 0;
}

static int print(
    const size_t step,
    const double time,
    const double dt
) {
  FILE * const stream = stdout;
  fprintf(
      stream,
      "step %10zu time % .4e dt % .4e\n",
      step, time, dt
  );
  return 0;
}

static int monitor_divergence(
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  const char file_name[] = ROOT_DIRECTORY "divergence.dat";
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hxxf = domain->hxxf;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  const double * const hzzf = domain->hzzf;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double div_max = 0.;
  double div_sum = 0.;
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double hx_xm = hxxf[i    ];
      const double hx_x0 = hxxc[i    ];
      const double hx_xp = hxxf[i + 1];
      const double hy_xm = hyxf[i    ];
      const double hy_x0 = hyxc[i    ];
      const double hy_xp = hyxf[i + 1];
      const double hz_zm = hzzf[k    ];
      const double hz_z0 = hzzc[k    ];
      const double hz_zp = hzzf[k + 1];
      const double ux_xm = ux[i    ][k    ];
      const double ux_xp = ux[i + 1][k    ];
      const double uz_zm = uz[i    ][k    ];
      const double uz_zp = uz[i    ][k + 1];
      const double div = 1. / jd(hx_x0, hy_x0, hz_z0) * (
          - jd(hx_xm, hy_xm, hz_z0) / hx_xm * ux_xm
          + jd(hx_xp, hy_xp, hz_z0) / hx_xp * ux_xp
          - jd(hx_x0, hy_x0, hz_zm) / hz_zm * uz_zm
          + jd(hx_x0, hy_x0, hz_zp) / hz_zp * uz_zp
      );
      div_max = fmax(div_max, fabs(div));
      div_sum = div_sum + div;
    }
  }
  return output(step, time, file_name, 2, (double []){div_max, div_sum});
}

static int monitor_max_velocity(
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  const char file_name[] = ROOT_DIRECTORY "max_velocity.dat";
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double ux_max = 0.;
  double uz_max = 0.;
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      ux_max = fmax(ux_max, fabs(ux[i][k]));
    }
  }
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      uz_max = fmax(uz_max, fabs(uz[i][k]));
    }
  }
  return output(step, time, file_name, 2, (double []){ux_max, uz_max});
}

static int monitor_interface_stats(
    const size_t step,
    const double time,
    const domain_t * const domain,
    const interface_field_t * const interface_field
) {
  const char file_name[] = ROOT_DIRECTORY "interface.dat";
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hyxc = domain->hyxc;
  const double * const hzzc = domain->hzzc;
  double ** const vof = interface_field->vof;
  double vof_sum = 0.;
  double vof_min = 0. + DBL_MAX;
  double vof_max = 0. - DBL_MAX;
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_vof = vof[i][k];
      vof_sum += l_vof * jd(hxxc[i], hyxc[i], hzzc[k]);
      vof_min = fmin(vof_min, l_vof);
      vof_max = fmax(vof_max, l_vof);
    }
  }
  return output(step, time, file_name, 3, (double []){vof_sum, vof_min, vof_max});
}

int monitor(
    const size_t step,
    const double time,
    const double dt,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field
) {
  if (0 != print(step, time, dt)) {
    LOGGER_FAILURE("failed to output metrics");
    goto abort;
  }
  if (0 != monitor_divergence(step, time, domain, flow_field)) {
    LOGGER_FAILURE("failed to check / output divergence");
    goto abort;
  }
  if (0 != monitor_max_velocity(step, time, domain, flow_field)) {
    LOGGER_FAILURE("failed to check / output maximum velocity");
    goto abort;
  }
  if (0 != monitor_interface_stats(step, time, domain, interface_field)) {
    LOGGER_FAILURE("failed to check / output interface statistics");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to monitor flow field");
  return 1;
}

