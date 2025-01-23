#include "./mass_flux.h"
#include "./compute_dux.h"

int compute_dux(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field,
    const double dt,
    double ** const force_x,
    double ** const dux
) {
  const size_t nx = domain->nx;
  const size_t nz = domain->nz;
  const double * const hxxc = domain->hxxc;
  const double * const hxxf = domain->hxxf;
  const double * const hyxc = domain->hyxc;
  const double * const hyxf = domain->hyxf;
  const double * const hzzc = domain->hzzc;
  double ** const ux = flow_field->ux;
  double ** const uz = flow_field->uz;
  double ** const p = flow_field->p;
  double ** const txx = flow_field->txx;
  double ** const tyy = flow_field->tyy;
  double ** const txz = flow_field->txz;
  double ** const vof_ux = interface_field->vof_ux;
  double ** const vof_uz = interface_field->vof_uz;
#pragma omp parallel for
  for (size_t i = ux_imin; i <= nx; i++) {
    for (size_t k = 1; k <= nz; k++) {
      const double l_jd = jd(hxxf[i], hyxf[i], hzzc[k]);
      const double mux_xm = + 0.5 * jdhx(hyxf[i - 1], hzzc[k    ]) * mass_flux(ux[i - 1][k    ], vof_ux[i - 1][k    ])
                            + 0.5 * jdhx(hyxf[i    ], hzzc[k    ]) * mass_flux(ux[i    ][k    ], vof_ux[i    ][k    ]);
      const double mux_xp = + 0.5 * jdhx(hyxf[i    ], hzzc[k    ]) * mass_flux(ux[i    ][k    ], vof_ux[i    ][k    ])
                            + 0.5 * jdhx(hyxf[i + 1], hzzc[k    ]) * mass_flux(ux[i + 1][k    ], vof_ux[i + 1][k    ]);
      const double muz_zm = + 0.5 * jdhz(hxxc[i - 1], hyxc[i - 1]) * mass_flux(uz[i - 1][k    ], vof_uz[i - 1][k    ])
                            + 0.5 * jdhz(hxxc[i    ], hyxc[i    ]) * mass_flux(uz[i    ][k    ], vof_uz[i    ][k    ]);
      const double muz_zp = + 0.5 * jdhz(hxxc[i - 1], hyxc[i - 1]) * mass_flux(uz[i - 1][k + 1], vof_uz[i - 1][k + 1])
                            + 0.5 * jdhz(hxxc[i    ], hyxc[i    ]) * mass_flux(uz[i    ][k + 1], vof_uz[i    ][k + 1]);
      const double ux_xm = + 0.5 * ux[i - 1][k    ]
                           + 0.5 * ux[i    ][k    ];
      const double ux_xp = + 0.5 * ux[i    ][k    ]
                           + 0.5 * ux[i + 1][k    ];
      const double ux_zm = + 0.5 * ux[i    ][k - 1]
                           + 0.5 * ux[i    ][k    ];
      const double ux_zp = + 0.5 * ux[i    ][k    ]
                           + 0.5 * ux[i    ][k + 1];
      const double advx = - 1. / l_jd * (
          - mux_xm * ux_xm
          + mux_xp * ux_xp
      );
      const double advz = - 1. / l_jd * (
          - muz_zm * ux_zm
          + muz_zp * ux_zp
      );
      const double difx = 1. / l_jd * (
          - jdhx(hyxc[i - 1], hzzc[k    ]) * txx[i - 1][k    ]
          + jdhx(hyxc[i    ], hzzc[k    ]) * txx[i    ][k    ]
      );
      const double difz = 1. / l_jd * (
          - jdhz(hxxf[i    ], hyxf[i    ]) * txz[i    ][k    ]
          + jdhz(hxxf[i    ], hyxf[i    ]) * txz[i    ][k + 1]
      );
      const double difa = - 1. / l_jd * (
          + 0.5 * (- jdhx(hyxf[i - 1], hzzc[k    ]) + jdhx(hyxf[i    ], hzzc[k    ])) * tyy[i - 1][k    ]
          + 0.5 * (- jdhx(hyxf[i    ], hzzc[k    ]) + jdhx(hyxf[i + 1], hzzc[k    ])) * tyy[i    ][k    ]
      );
      const double pre = - 1. / hxxf[i] * (
          - p[i - 1][k    ]
          + p[i    ][k    ]
      );
      const double surf = force_x[i][k];
      dux[i][k] = dt * (
          + advx + advz
          + difx + difz + difa
          + pre
          + surf
      );
    }
  }
  return 0;
}

