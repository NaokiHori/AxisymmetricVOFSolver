#include <stdlib.h>
#include <math.h>
#include "logger.h"
#include "memory.h"
#include "visualize.h"
#include "./visualize/contouring.h"
#include "./visualize/drawer.h"

int visualize(
    const char file_name[],
    const size_t nx,
    const size_t nz,
    const double * const x,
    const double * const z,
    const double * const * const vof
) {
  // compute iso-lines
  contouring_contour_t * contouring_contour = NULL;
  if (0 != contouring_exec(0.5, (size_t []){nz + 2, nx + 2}, (const double * const []){z, x}, vof, &contouring_contour)) {
    LOGGER_FAILURE("failed to extract iso-lines");
    return 1;
  }
  // convert data type (from singly-linked lists to normal arrays)
  //   for later convenience
  size_t ncontours = 0;
  for (const contouring_contour_t * c = contouring_contour; c; c = c->next) {
    ncontours += 1;
  }
  contour_t * const contours = memory_alloc(ncontours, sizeof(contour_t));
  {
    size_t contour_id = 0;
    for (const contouring_contour_t * c = contouring_contour; c; c = c->next) {
      const size_t npoints = c->npoints;
      contour_t * const contour = contours + contour_id;
      contour->is_closed = c->is_closed;
      contour->npoints = npoints;
      contour->points = memory_alloc(npoints, sizeof(point_t));
      size_t point_id = 0;
      for (const contouring_point_t * p = c->point; p; p = p->next) {
        point_t * const point = contour->points + point_id;
        point->x = p->x;
        point->y = p->y;
        point->nx = nan("");
        point->ny = nan("");
        point_id += 1;
      }
      contour_id += 1;
    }
  }
  if (0 != contouring_cleanup(&contouring_contour)) {
    return 1;
  }
  // compute vertex normal
  for (size_t contour_id = 0; contour_id < ncontours; contour_id++) {
    contour_t * const contour = contours + contour_id;
    const size_t npoints = contour->npoints;
    point_t * const points = contour->points;
    if (contour->is_closed) {
      for (size_t point_id = 0; point_id < npoints; point_id++) {
        point_t * const point_m = points + (point_id - 1 + npoints) % npoints;
        point_t * const point_0 = points + (point_id     + npoints) % npoints;
        point_t * const point_p = points + (point_id + 1 + npoints) % npoints;
        const double distance_m = sqrt(
            + pow(point_m->x - point_0->x, 2)
            + pow(point_m->y - point_0->y, 2)
        );
        const double distance_p = sqrt(
            + pow(point_0->x - point_p->x, 2)
            + pow(point_0->y - point_p->y, 2)
        );
        const double dxdt = 1. / (distance_m + distance_p) * (
            + distance_m / distance_p * (point_p->x - point_0->x)
            - distance_p / distance_m * (point_m->x - point_0->x)
        );
        const double dydt = 1. / (distance_m + distance_p) * (
            + distance_m / distance_p * (point_p->y - point_0->y)
            - distance_p / distance_m * (point_m->y - point_0->y)
        );
        const double norm = sqrt(pow(dxdt, 2.) + pow(dydt, 2.));
        point_0->nx = + dydt / norm;
        point_0->ny = - dxdt / norm;
      }
    } else {
      for (size_t point_id = 0; point_id < npoints; point_id++) {
        if (0 == point_id) {
          point_t * const point_0 = points + 0;
          point_t * const point_p = points + 1;
          const double distance = sqrt(
              + pow(point_p->x - point_0->x, 2)
              + pow(point_p->y - point_0->y, 2)
          );
          const double dxdt = (point_p->x - point_0->x) / distance;
          const double dydt = (point_p->y - point_0->y) / distance;
          const double norm = sqrt(pow(dxdt, 2.) + pow(dydt, 2.));
          point_0->nx = + dydt / norm;
          point_0->ny = - dxdt / norm;
        } else if (npoints - 1 == point_id) {
          point_t * const point_m = points + npoints - 2;
          point_t * const point_0 = points + npoints - 1;
          const double distance = sqrt(
              + pow(point_m->x - point_0->x, 2)
              + pow(point_m->y - point_0->y, 2)
          );
          const double dxdt = (point_0->x - point_m->x) / distance;
          const double dydt = (point_0->y - point_m->y) / distance;
          const double norm = sqrt(pow(dxdt, 2.) + pow(dydt, 2.));
          point_0->nx = + dydt / norm;
          point_0->ny = - dxdt / norm;
        } else {
          point_t * const point_m = points + point_id - 1;
          point_t * const point_0 = points + point_id    ;
          point_t * const point_p = points + point_id + 1;
          const double distance_m = sqrt(
              + pow(point_m->x - point_0->x, 2)
              + pow(point_m->y - point_0->y, 2)
          );
          const double distance_p = sqrt(
              + pow(point_0->x - point_p->x, 2)
              + pow(point_0->y - point_p->y, 2)
          );
          const double dxdt = 1. / (distance_m + distance_p) * (
              + distance_m / distance_p * (point_p->x - point_0->x)
              - distance_p / distance_m * (point_m->x - point_0->x)
          );
          const double dydt = 1. / (distance_m + distance_p) * (
              + distance_m / distance_p * (point_p->y - point_0->y)
              - distance_p / distance_m * (point_m->y - point_0->y)
          );
          const double norm = sqrt(pow(dxdt, 2.) + pow(dydt, 2.));
          point_0->nx = + dydt / norm;
          point_0->ny = - dxdt / norm;
        }
      }
    }
  }
  // render three-dimensional scene
  if (0 != draw(file_name, ncontours, contours)) {
    LOGGER_FAILURE("failed to draw image");
    return 1;
  }
  // clean-up
  for (size_t contour_id = 0; contour_id < ncontours; contour_id++) {
    const contour_t * const contour = contours + contour_id;
    memory_free(contour->points);
  }
  memory_free(contours);
  return 0;
}

