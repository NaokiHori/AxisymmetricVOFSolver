#include "param.h"

double mass_flux(
    const double velocity,
    const double vof
) {
  return velocity * (1. + (DENSITY_RATIO - 1.) * vof);
}
