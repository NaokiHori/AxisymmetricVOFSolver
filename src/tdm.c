#include <stdio.h>
#include <stdbool.h>
#include <float.h> // DBL_EPSILON
#include "memory.h"
#include "tdm.h"

struct tdm_internal_t {
  // auxiliary buffers
  double * v;
  double * w;
};

int tdm_init_plan(
    const size_t nitems,
    const size_t repeat_for,
    const bool is_periodic,
    tdm_plan_t ** tdm_plan
) {
  const size_t minimum_nitems = 3;
  if (nitems < minimum_nitems) {
    fprintf(stderr, "size %zu is too small, give larger than %zu\n", nitems, minimum_nitems);
    return 1;
  }
  *tdm_plan = memory_alloc(1, sizeof(tdm_plan_t));
  (*tdm_plan)->internal = memory_alloc(1, sizeof(tdm_internal_t));
  (*tdm_plan)->nitems = nitems;
  (*tdm_plan)->repeat_for = repeat_for;
  (*tdm_plan)->is_periodic = is_periodic;
  (*tdm_plan)->internal->v = memory_alloc(nitems * repeat_for, sizeof(double));
  (*tdm_plan)->internal->w = memory_alloc(nitems * repeat_for, sizeof(double));
  return 0;
}

static double myfabs(
    const double v
) {
  return v < 0. ? - v : v;
}

int tdm_solve(
    tdm_plan_t * const tdm_plan,
    const double * const l,
    const double * const c,
    const double * const u,
    const double * const c_offsets,
    double * const qs
) {
  if (NULL == tdm_plan) {
    return 1;
  }
  const size_t nitems = tdm_plan->nitems;
  const size_t repeat_for = tdm_plan->repeat_for;
  const bool is_periodic = tdm_plan->is_periodic;
#pragma omp parallel for
  for (size_t j = 0; j < repeat_for; j++) {
    const double c_offset = c_offsets[j];
    double * const v = tdm_plan->internal->v + j * nitems;
    double * const w = tdm_plan->internal->w + j * nitems;
    double * const q = qs + j * nitems;
    if (is_periodic) {
      // consider a perturbed system as well
      for (size_t i = 0; i < nitems - 1; i++) {
        w[i]
          = i ==          0 ? - 1. * l[i]
          : i == nitems - 2 ? - 1. * u[i]
          : 0.;
      }
      // divide the first row by center-diagonal term
      v[0] = u[0] / (c[0] + c_offset);
      q[0] = q[0] / (c[0] + c_offset);
      w[0] = w[0] / (c[0] + c_offset);
      // forward sweep
      for (size_t i = 1; i < nitems - 1; i++) {
        // assume positive-definite system
        //   to skip zero-division checks
        const double val = 1. / (c[i] + c_offset - l[i] * v[i-1]);
        v[i] = val * u[i];
        q[i] = val * (q[i] - l[i] * q[i-1]);
        w[i] = val * (w[i] - l[i] * w[i-1]);
      }
      // backward substitution
      for (size_t i = nitems - 3; ; i--) {
        q[i] -= v[i] * q[i+1];
        w[i] -= v[i] * w[i+1];
        if (0 == i) {
          break;
        }
      }
      // couple two systems to find the answer
      const double num = q[nitems - 1]            - u[nitems - 1] * q[0] - l[nitems - 1] * q[nitems - 2];
      const double den = c[nitems - 1] + c_offset + u[nitems - 1] * w[0] + l[nitems - 1] * w[nitems - 2];
      q[nitems - 1] = myfabs(den) < DBL_EPSILON ? 0. : num / den;
      for (size_t i = 0; i < nitems - 1; i++) {
        q[i] = q[i] + q[nitems - 1] * w[i];
      }
    } else {
      // divide the first row by center-diagonal term
      v[0] = u[0] / (c[0] + c_offset);
      q[0] = q[0] / (c[0] + c_offset);
      // forward sweep
      for (size_t i = 1; i < nitems - 1; i++) {
        // assume positive-definite system
        //   to skip zero-division checks
        const double val = 1. / (c[i] + c_offset - l[i] * v[i - 1]);
        v[i] = val * u[i];
        q[i] = val * (q[i] - l[i] * q[i - 1]);
      }
      // last row, do the same thing but consider singularity (degeneracy)
      const double val = c[nitems - 1] + c_offset - l[nitems - 1] * v[nitems - 2];
      if (DBL_EPSILON < myfabs(val)) {
        q[nitems - 1] = 1. / val * (q[nitems - 1] - l[nitems - 1] * q[nitems - 2]);
      } else {
        // singular
        q[nitems - 1] = 0.;
      }
      // backward substitution
      for (size_t i = nitems - 2; ; i--) {
        q[i] -= v[i] * q[i + 1];
        if (0 == i) {
          break;
        }
      }
    }
  }
  return 0;
}

int tdm_destroy_plan(
    tdm_plan_t ** tdm_plan
) {
  memory_free((*tdm_plan)->internal->v);
  memory_free((*tdm_plan)->internal->w);
  memory_free((*tdm_plan)->internal);
  memory_free(*tdm_plan);
  *tdm_plan = NULL;
  return 0;
}

