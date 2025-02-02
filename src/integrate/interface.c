#include <math.h>
#include "param.h"
#include "domain.h"
#include "impose_bc.h"
#include "exchange_halo.h"
#include "./interface.h"

#define N_GAUSS 2

static const double gauss_ws[N_GAUSS] = {
  0.5,
  0.5,
};

static const double gauss_ps[N_GAUSS] = {
  - 0.2886751345948129,
  + 0.2886751345948129,
};

// a small number to avoid zero-division
static const double small = 1.e-16;

static const double vof_min = 1.e-8;

static double indicator(
    const double vof_beta,
    const double nx,
    const double nz,
    const double sdf,
    const double x,
    const double z
) {
  const double surface_function = nx * x + nz * z + sdf;
  return 0.5 * (1. + tanh(vof_beta * surface_function));
}

static double vof_to_sdf(
    const double vof_beta,
    const double vof
) {
  const double clamped_vof = fmin(1. - small, fmax(0. + small, vof));
  return - 0.5 / vof_beta * log(1. / clamped_vof - 1.);
}

static double compute_norm(
    const double x,
    const double z
) {
  return fmax(
      small,
      sqrt(
        + pow(x, 2.)
        + pow(z, 2.)
      )
  );
}

static int compute_metrics(
    const domain_t * const domain,
    interface_field_t * const interface_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hxxf = domain->hxxf;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  const double * const hzzf = domain->hzzf;
  const double * const * const vof = (const double * const * const)interface_field->vof;
  double * const * const corner_nx = interface_field->corner_nx;
  double * const * const corner_nz = interface_field->corner_nz;
  double * const * const center_nx = interface_field->center_nx;
  double * const * const center_nz = interface_field->center_nz;
  double * const * const curv = interface_field->curv;
#pragma omp parallel for
  for (size_t i = 1; i <= nx + 1; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      double * const l_corner_nx = &corner_nx[i][k];
      double * const l_corner_nz = &corner_nz[i][k];
      *l_corner_nx = 0.5 * (
          - vof[i - 1][k - 1] + vof[i    ][k - 1]
          - vof[i - 1][k    ] + vof[i    ][k    ]
      ) / hxxf[i    ];
      *l_corner_nz = 0.5 * (
          - vof[i - 1][k - 1] - vof[i    ][k - 1]
          + vof[i - 1][k    ] + vof[i    ][k    ]
      ) / hzzf[k    ];
      const double norm = compute_norm(*l_corner_nx, *l_corner_nz);
      *l_corner_nx /= norm;
      *l_corner_nz /= norm;
    }
  }
  if (Z_PERIODIC) {
    exchange_halo_z(domain, corner_nx);
    exchange_halo_z(domain, corner_nz);
  } else {
    for (size_t i = 0; i <= nx + 1; i++) {
      corner_nx[i][     0] = 0.; // unused
      corner_nx[i][     1] = 0.;
      corner_nx[i][nz + 1] = 0.;
      corner_nz[i][     0] = 0.; // unused
      corner_nz[i][     1] = 0.;
      corner_nz[i][nz + 1] = 0.;
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double hx_x0 = hxxc[i    ];
      const double hy_xm = hyxf[i    ];
      const double hy_x0 = hyxc[i    ];
      const double hy_xp = hyxf[i + 1];
      const double hz_z0 = hzzc[k    ];
      const double nx_xm =
        + 0.5 * corner_nx[i    ][k    ]
        + 0.5 * corner_nx[i    ][k + 1];
      const double nx_xp =
        + 0.5 * corner_nx[i + 1][k    ]
        + 0.5 * corner_nx[i + 1][k + 1];
      const double nz_zm =
        + 0.5 * corner_nz[i    ][k    ]
        + 0.5 * corner_nz[i + 1][k    ];
      const double nz_zp =
        + 0.5 * corner_nz[i    ][k + 1]
        + 0.5 * corner_nz[i + 1][k + 1];
      curv[i][k] = - 1. / jd(hx_x0, hy_x0, hz_z0) * (
          - jdhx(hy_xm, hz_z0) * nx_xm
          + jdhx(hy_xp, hz_z0) * nx_xp
          - jdhz(hx_x0, hy_x0) * nz_zm
          + jdhz(hx_x0, hy_x0) * nz_zp
      );
      center_nx[i][k] =
        + 0.25 * corner_nx[i    ][k    ]
        + 0.25 * corner_nx[i + 1][k    ]
        + 0.25 * corner_nx[i    ][k + 1]
        + 0.25 * corner_nx[i + 1][k + 1];
      center_nz[i][k] =
        + 0.25 * corner_nz[i    ][k    ]
        + 0.25 * corner_nz[i + 1][k    ]
        + 0.25 * corner_nz[i    ][k + 1]
        + 0.25 * corner_nz[i + 1][k + 1];
      const double norm = compute_norm(center_nx[i][k], center_nz[i][k]);
      center_nx[i][k] /= norm;
      center_nz[i][k] /= norm;
    }
  }
  if (Z_PERIODIC) {
    exchange_halo_z(domain, curv);
    exchange_halo_z(domain, center_nx);
    exchange_halo_z(domain, center_nz);
  }
  return 0;
}

static int compute_force(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    interface_field_t * const interface_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxf = domain->hxxf;
  const double * const hzzf = domain->hzzf;
  double ** const density = flow_field->density[0];
  double ** const vof = interface_field->vof;
  double ** const curv = interface_field->curv;
  double ** const force_x = interface_field->force_x;
  double ** const force_z = interface_field->force_z;
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[i - 1][k    ]
        + 0.5 * density[i    ][k    ];
      const double l_curv =
        + 0.5 * curv[i - 1][k    ]
        + 0.5 * curv[i    ][k    ];
      const double l_dvofdx = (
          - vof[i - 1][k    ]
          + vof[i    ][k    ]
      ) / hxxf[i    ];
      force_x[i][k] = 2. * l_density / (1. + DENSITY_RATIO) / We * l_curv * l_dvofdx;
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      const double l_density =
        + 0.5 * density[i    ][k - 1]
        + 0.5 * density[i    ][k    ];
      const double l_curv =
        + 0.5 * curv[i    ][k - 1]
        + 0.5 * curv[i    ][k    ];
      const double l_dvofdz = (
          - vof[i    ][k - 1]
          + vof[i    ][k    ]
      ) / hzzf[k    ];
      force_z[i][k] = 2. * l_density / (1. + DENSITY_RATIO) / We * l_curv * l_dvofdz;
    }
  }
  return 0;
}

static int compute_face_vof(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    interface_field_t * const interface_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double vof_beta = interface_field->vof_beta;
  double ** const vof = interface_field->vof;
  double ** const center_nx = interface_field->center_nx;
  double ** const center_nz = interface_field->center_nz;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const vof_uz = interface_field->vof_uz;
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_vel = ux[i][k];
      const size_t ii = l_vel < 0. ? i : i - 1;
      const double x = l_vel < 0. ? -0.5 : 0.5;
      const double l_vof = vof[ii][k];
      double * const l_vof_ux = &vof_ux[i][k];
      if (l_vof < vof_min || 1. - vof_min < l_vof) {
        *l_vof_ux = l_vof;
      } else {
        const double l_center_nx = center_nx[ii][k];
        const double l_center_nz = center_nz[ii][k];
        const double l_sdf = vof_to_sdf(vof_beta, l_vof);
        *l_vof_ux = 0.;
        for (size_t kk = 0; kk < N_GAUSS; kk++) {
          const double w = gauss_ws[kk];
          const double z = gauss_ps[kk];
          *l_vof_ux += w * indicator(vof_beta, l_center_nx, l_center_nz, l_sdf, x, z);
        }
      }
    }
  }
  for (size_t k = 0; k <= nz + 1; k++) {
    vof_ux[     0][k] = nan(""); // unused
    vof_ux[     1][k] = vof[     0][k];
    vof_ux[nx + 1][k] = vof[nx + 1][k];
  }
  if (Z_PERIODIC) {
    exchange_halo_z(domain, vof_ux);
  } else {
    for (size_t i = 1; i <= nx; i++) {
      vof_ux[i][     0] = vof[i][     0];
      vof_ux[i][nz + 1] = vof[i][nz + 1];
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      const double l_vel = uz[i][k];
      const size_t kk = l_vel < 0. ? k : k - 1;
      const double z = l_vel < 0. ? -0.5 : 0.5;
      const double l_vof = vof[i][kk];
      double * const l_vof_uz = &vof_uz[i][k];
      if (l_vof < vof_min || 1. - vof_min < l_vof) {
        *l_vof_uz = l_vof;
      } else {
        const double l_center_nx = center_nx[i][kk];
        const double l_center_nz = center_nz[i][kk];
        const double l_sdf = vof_to_sdf(vof_beta, l_vof);
        *l_vof_uz = 0.;
        for (size_t ii = 0; ii < N_GAUSS; ii++) {
          const double w = gauss_ws[ii];
          const double x = gauss_ps[ii];
          *l_vof_uz += w * indicator(vof_beta, l_center_nx, l_center_nz, l_sdf, x, z);
        }
      }
    }
  }
  for (size_t k = 0; k <= nz + 1; k++) {
    vof_uz[     0][k] = vof[     0][k];
    vof_uz[nx + 1][k] = vof[nx + 1][k];
  }
  if (Z_PERIODIC) {
    exchange_halo_z(domain, vof_uz);
  } else {
    for (size_t i = 1; i <= nx; i++) {
      vof_uz[i][     0] = nan(""); // unused
      vof_uz[i][     1] = vof[i][     0];
      vof_uz[i][nz + 1] = vof[i][nz + 1];
    }
  }
  return 0;
}

static int update_vof(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    interface_field_t * const interface_field,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const vof_uz = interface_field->vof_uz;
  double ** const vof = interface_field->vof;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double hx_x0 = hxxc[i    ];
      const double hy_xm = hyxf[i    ];
      const double hy_x0 = hyxc[i    ];
      const double hy_xp = hyxf[i + 1];
      const double hz_z0 = hzzc[k    ];
      const double flux_xm = ux[i    ][k    ] * vof_ux[i    ][k    ];
      const double flux_xp = ux[i + 1][k    ] * vof_ux[i + 1][k    ];
      const double flux_zm = uz[i    ][k    ] * vof_uz[i    ][k    ];
      const double flux_zp = uz[i    ][k + 1] * vof_uz[i    ][k + 1];
      const double dvof = 1. / jd(hx_x0, hy_x0, hz_z0) * (
          + jdhx(hy_xm, hz_z0) * flux_xm
          - jdhx(hy_xp, hz_z0) * flux_xp
          + jdhz(hx_x0, hy_x0) * flux_zm
          - jdhz(hx_x0, hy_x0) * flux_zp
      );
      vof[i][k] += dvof * dt;
    }
  }
  impose_bc_vof_x(domain, vof);
  if (Z_PERIODIC) {
    exchange_halo_z(domain, vof);
  } else {
    impose_bc_vof_z(domain, vof);
  }
  return 0;
}

int update_interface(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    interface_field_t * const interface_field,
    const double dt
) {
  compute_metrics(domain, interface_field);
  compute_force(domain, flow_field, interface_field);
  compute_face_vof(domain, flow_field, interface_field);
  update_vof(domain, flow_field, interface_field, dt);
  return 0;
}

