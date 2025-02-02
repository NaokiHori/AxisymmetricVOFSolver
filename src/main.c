#include <stdio.h>
#include <stddef.h> // size_t
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "interface_field.h"
#include "integrate.h"
#include "monitor.h"
#include "save.h"
#include "./visualize.h"

typedef struct {
  double monitor;
  double save;
  double visualize;
} schedule_t;

int main(
    void
) {
  domain_t domain = {};
  flow_field_t flow_field = {};
  flow_solver_t flow_solver = {};
  interface_field_t interface_field = {};
  if (0 != domain_init(&domain)) {
    return 1;
  }
  if (0 != flow_field_init(&domain, &flow_field)) {
    return 1;
  }
  if (0 != flow_solver_init(&domain, &flow_solver)) {
    return 1;
  }
  if (0 != interface_field_init(&domain, &interface_field)) {
    return 1;
  }
  const double time_max = 1.5e+1;
  const schedule_t rate = {
    .monitor = 1.e-1,
    .save = 1.e+0,
    .visualize = 5.e-2,
  };
  schedule_t next = rate;
  for (double time = 0.; time < time_max; ) {
    static size_t step = 0;
    double dt = 0.;
    if (0 != integrate(&domain, &flow_field, &flow_solver, &interface_field, &dt)) {
      break;
    }
    step += 1;
    time += dt;
    if (next.monitor < time) {
      monitor(step, time, dt, &domain, &flow_field, &interface_field);
      next.monitor += rate.monitor;
    }
    if (next.save < time) {
      static size_t id = 0;
      save(id, step, time, &domain, &flow_field, &interface_field);
      id += 1;
      next.save += rate.save;
    }
    if (next.visualize < time) {
      static size_t id = 0;
      char file_name[256] = {'\0'};
      snprintf(file_name, sizeof(file_name) - 1, "output/image/image%010zu.pgm", id);
      visualize(file_name, domain.nx, domain.nz, domain.xc, domain.zc, (const double * const *)interface_field.vof);
      id += 1;
      next.visualize += rate.visualize;
    }
  }
  if (0 != domain_finalize(&domain)) {
    return 1;
  }
  if (0 != flow_field_finalize(&flow_field)) {
    return 1;
  }
  if (0 != flow_solver_finalize(&flow_solver)) {
    return 1;
  }
  if (0 != interface_field_finalize(&interface_field)) {
    return 1;
  }
  return 0;
}

