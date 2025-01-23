#include "param.h"
#include "./mass_flux.h"
#include "./compute_duz.h"

int compute_duz(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field,
    const double dt,
    double ** const force_z,
    double ** const duz
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  const double * const hzzf = domain->hzzf;
  double ** const density = flow_field->density[0];
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double ** const p = flow_field->p;
  double ** const txz = flow_field->txz;
  double ** const tzz = flow_field->tzz;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const vof_uz = interface_field->vof_uz;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t k = uz_kmin; k <= nz; k++) {
      const double l_jd = jd(hxxc[i], hyxc[i], hzzf[k]);
      const double mux_xm = + 0.5 * jdhx(hyxf[i    ], hzzc[k - 1]) * mass_flux(ux[i    ][k - 1], vof_ux[i    ][k - 1])
                            + 0.5 * jdhx(hyxf[i    ], hzzc[k    ]) * mass_flux(ux[i    ][k    ], vof_ux[i    ][k    ]);
      const double mux_xp = + 0.5 * jdhx(hyxf[i + 1], hzzc[k - 1]) * mass_flux(ux[i + 1][k - 1], vof_ux[i + 1][k - 1])
                            + 0.5 * jdhx(hyxf[i + 1], hzzc[k    ]) * mass_flux(ux[i + 1][k    ], vof_ux[i + 1][k    ]);
      const double muz_zm = + 0.5 * jdhz(hxxc[i    ], hyxc[i    ]) * mass_flux(uz[i    ][k - 1], vof_uz[i    ][k - 1])
                            + 0.5 * jdhz(hxxc[i    ], hyxc[i    ]) * mass_flux(uz[i    ][k    ], vof_uz[i    ][k    ]);
      const double muz_zp = + 0.5 * jdhz(hxxc[i    ], hyxc[i    ]) * mass_flux(uz[i    ][k    ], vof_uz[i    ][k    ])
                            + 0.5 * jdhz(hxxc[i    ], hyxc[i    ]) * mass_flux(uz[i    ][k + 1], vof_uz[i    ][k + 1]);
      const double uz_xm = + 0.5 * uz[i - 1][k    ]
                           + 0.5 * uz[i    ][k    ];
      const double uz_xp = + 0.5 * uz[i    ][k    ]
                           + 0.5 * uz[i + 1][k    ];
      const double uz_zm = + 0.5 * uz[i    ][k - 1]
                           + 0.5 * uz[i    ][k    ];
      const double uz_zp = + 0.5 * uz[i    ][k    ]
                           + 0.5 * uz[i    ][k + 1];
      const double advx = - 1. / l_jd * (
          - mux_xm * uz_xm
          + mux_xp * uz_xp
      );
      const double advz = - 1. / l_jd * (
          - muz_zm * uz_zm
          + muz_zp * uz_zp
      );
      const double difx = 1. / l_jd * (
          - jdhx(hyxf[i    ], hzzf[k    ]) * txz[i    ][k    ]
          + jdhx(hyxf[i + 1], hzzf[k    ]) * txz[i + 1][k    ]
      );
      const double difz = 1. / l_jd * (
          - jdhz(hxxc[i    ], hyxc[i    ]) * tzz[i    ][k - 1]
          + jdhz(hxxc[i    ], hyxc[i    ]) * tzz[i    ][k    ]
      );
      const double pre = - 1. / hzzf[k] * (
          - p[i    ][k - 1]
          + p[i    ][k    ]
      );
      const double surf = force_z[i][k];
      const double gravity = - 1. / Fr / Fr * (
          + 0.5 * density[i    ][k - 1]
          + 0.5 * density[i    ][k    ]
      );
      duz[i][k] = dt * (
          + advx + advz
          + difx + difz
          + pre
          + surf
          + gravity
      );
    }
  }
  return 0;
}

