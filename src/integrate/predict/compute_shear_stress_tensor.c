#include <math.h>
#include "logger.h"
#include "param.h"
#include "exchange_halo.h"
#include "compute_shear_stress_tensor.h"

static double convert(
    double vof
) {
  vof = vof < 0. ? 0. : vof;
  vof = 1. < vof ? 1. : vof;
  return 1. / Re * (1. + (VISCOSITY_RATIO - 1.) * vof);
}

static int compute_txx(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  double ** const ux = flow_field->ux;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const txx = flow_field->txx;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_vis =
        + 0.5 * convert(vof_ux[i    ][k    ])
        + 0.5 * convert(vof_ux[i + 1][k    ]);
      const double lxx = 1. / hxxc[i    ] * (
          - ux[i    ][k    ]
          + ux[i + 1][k    ]
      );
      txx[i][k] = 2. * l_vis * lxx;
    }
  }
  for (size_t k = 0; k <= nz + 1; k++) {
    txx[     0][k] = nan("");
    txx[nx + 1][k] = nan("");
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, txx)) {
      LOGGER_FAILURE("failed to exchange halo in z (txx)");
      goto abort;
    }
  } else {
    for (size_t i = 0; i <= nx + 1; i++) {
      txx[i][     0] = nan("");
      txx[i][nz + 1] = nan("");
    }
  }
  return 0;
abort:
  return 1;
}

static int compute_tyy(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hxxf = domain->hxxf;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  double ** const ux = flow_field->ux;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const tyy = flow_field->tyy;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_vis =
        + 0.5 * convert(vof_ux[i    ][k    ])
        + 0.5 * convert(vof_ux[i + 1][k    ]);
      const double lyy = 1. / jd(hxxc[i    ], hyxc[i    ], hzzc[k    ]) * (
          - jd(hxxf[i    ], hyxf[i    ], hzzc[k    ]) / hxxf[i    ]
          + jd(hxxf[i + 1], hyxf[i + 1], hzzc[k    ]) / hxxf[i + 1]
      ) * (
          + 0.5 * ux[i    ][k    ]
          + 0.5 * ux[i + 1][k    ]
      );
      tyy[i][k] = 2. * l_vis * lyy;
    }
  }
  for (size_t k = 0; k <= nz + 1; k++) {
    tyy[     0][k] = nan("");
    tyy[nx + 1][k] = nan("");
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, tyy)) {
      LOGGER_FAILURE("failed to exchange halo in z (tyy)");
      goto abort;
    }
  } else {
    for (size_t i = 0; i <= nx + 1; i++) {
      tyy[i][     0] = nan("");
      tyy[i][nz + 1] = nan("");
    }
  }
  return 0;
abort:
  return 1;
}

static int compute_txz(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxf = domain->hxxf;
  const double * const hzzf = domain->hzzf;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const vof_uz = interface_field->vof_uz;
  double ** const txz = flow_field->txz;
#pragma omp parallel for
  for (size_t i = 1; i <= nx + 1; i++) {
    for (size_t k = 1; k <= nz + 1; k++) {
      const double l_vis_ux =
        + 0.5 * convert(vof_ux[i    ][k - 1])
        + 0.5 * convert(vof_ux[i    ][k    ]);
      const double l_vis_uz =
        + 0.5 * convert(vof_uz[i - 1][k    ])
        + 0.5 * convert(vof_uz[i    ][k    ]);
      const double lxz = 1. / hzzf[k    ] * (
          - ux[i    ][k - 1]
          + ux[i    ][k    ]
      );
      const double lzx = 1. / hxxf[i    ] * (
          - uz[i - 1][k    ]
          + uz[i    ][k    ]
      );
      txz[i][k] = l_vis_ux * lxz + l_vis_uz * lzx;
    }
  }
  for (size_t k = 0; k <= nz + 1; k++) {
    txz[0][k] = nan("");
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, txz)) {
      LOGGER_FAILURE("failed to exchange halo in z (txz)");
      goto abort;
    }
  } else {
    for (size_t i = 0; i <= nx + 1; i++) {
      txz[i][0] = nan("");
    }
  }
  return 0;
abort:
  return 1;
}

static int compute_tzz(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hzzc = domain->hzzc;
  double ** const uz = flow_field->uz;
  double ** const vof_uz = interface_field->vof_uz;
  double ** const tzz = flow_field->tzz;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_vis =
        + 0.5 * convert(vof_uz[i    ][k    ])
        + 0.5 * convert(vof_uz[i    ][k + 1]);
      const double lzz = 1. / hzzc[k    ] * (
          - uz[i    ][k    ]
          + uz[i    ][k + 1]
      );
      tzz[i][k] = 2. * l_vis * lzz;
    }
  }
  for (size_t k = 0; k <= nz + 1; k++) {
    tzz[     0][k] = nan("");
    tzz[nx + 1][k] = nan("");
  }
  if (Z_PERIODIC) {
    if (0 != exchange_halo_z(domain, tzz)) {
      LOGGER_FAILURE("failed to exchange halo in z (tzz)");
      goto abort;
    }
  } else {
    for (size_t i = 0; i <= nx + 1; i++) {
      tzz[i][     0] = nan("");
      tzz[i][nz + 1] = nan("");
    }
  }
  return 0;
abort:
  return 1;
}

int compute_shear_stress_tensor(
    const domain_t * const domain,
    const interface_field_t * const interface_field,
    flow_field_t * const flow_field
) {
  if (0 != compute_txx(domain, interface_field, flow_field)) {
    LOGGER_FAILURE("failed to compute txx");
    goto abort;
  }
  if (0 != compute_tyy(domain, interface_field, flow_field)) {
    LOGGER_FAILURE("failed to compute tyy");
    goto abort;
  }
  if (0 != compute_txz(domain, interface_field, flow_field)) {
    LOGGER_FAILURE("failed to compute txz");
    goto abort;
  }
  if (0 != compute_tzz(domain, interface_field, flow_field)) {
    LOGGER_FAILURE("failed to compute tzz");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

