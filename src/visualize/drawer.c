#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include "./drawer.h"
#include "./drawer/types.h"
#include "./drawer/vector.h"
#include "./drawer/render.h"
#include "./drawer/output.h"

static const double pi = 3.1415926535897932384626;

static vector_t rotate_vector(
    const vector_t normal,
    const double theta,
    const vector_t target
) {
  const double nx = normal.x;
  const double ny = normal.y;
  const double nz = normal.z;
  const double c = cos(theta);
  const double s = sin(theta);
  const double matrix[3][3] = {
    {
      nx * nx * (1. - c) + c,
      nx * ny * (1. - c) - nz * s,
      nx * nz * (1. - c) + ny * s,
    },
    {
      ny * nx * (1. - c) + nz * s,
      ny * ny * (1. - c) + c,
      ny * nz * (1. - c) - nx * s,
    },
    {
      nz * nx * (1. - c) - ny * s,
      nz * ny * (1. - c) + nx * s,
      nz * nz * (1. - c) + c,
    },
  };
  return (vector_t){
    .x = matrix[0][0] * target.x + matrix[0][1] * target.y + matrix[0][2] * target.z,
    .y = matrix[1][0] * target.x + matrix[1][1] * target.y + matrix[1][2] * target.z,
    .z = matrix[2][0] * target.x + matrix[2][1] * target.y + matrix[2][2] * target.z,
  };
}

int draw(
    const char file_name[],
    const size_t ncontours,
    const contour_t * const contours
) {
  const vector_t elevation_vector = {
    .x = 1.,
    .y = 0.,
    .z = 0.,
  };
  const double elevation_angle = - pi / 12.;
  const camera_t camera = {
    .position = rotate_vector(
        elevation_vector,
        elevation_angle,
        (vector_t){.x = 0., .y = -16., .z = 0.}
    ),
    .look_at = {.x = 0., .y = 0., .z = 0.},
  };
  const vector_t light = normalize(mul(-1., camera.position));
  const size_t screen_resolution[2] = {768, 1536};
  const vector_t screen_center = rotate_vector(
      elevation_vector,
      elevation_angle,
      (vector_t){.x = 0., .y = -10., .z = 0.}
  );
  const vector_t screen_local_x = rotate_vector(
      elevation_vector,
      elevation_angle,
      (vector_t){.x = 2., .y = 0., .z = 0.}
  );
  const vector_t screen_local_y = rotate_vector(
      elevation_vector,
      elevation_angle,
      (vector_t){.x = 0., .y = 0., .z = 4.}
  );
  const screen_t screen = {
    .center = screen_center,
    .width = screen_resolution[0],
    .height = screen_resolution[1],
    .local_x = screen_local_x,
    .local_y = screen_local_y,
    .normal = normalize(
        outer_product(screen_local_x, screen_local_y)
    ),
  };
  pixel_t * const canvas = malloc(screen.width * screen.height * sizeof(pixel_t));
  for (size_t n = 0; n < screen.width * screen.height; n++) {
    pixel_t * const pixel = canvas + n;
    pixel->color = (color_t){
      .grayscale = 0x00,
    };
    pixel->depth = -1. * DBL_MAX;
  }
  for (size_t contour_id = 0; contour_id < ncontours; contour_id++) {
    const contour_t * const c = contours + contour_id;
    const size_t npoints = c->npoints;
    point_t * const points = c->points;
    const size_t n_azimuth = 128;
    for (size_t m = 0; m < n_azimuth; m++) {
      const double theta0 = 2. * pi * (m + 0) / n_azimuth;
      const double theta1 = 2. * pi * (m + 1) / n_azimuth;
      const double cos_theta0 = cos(theta0);
      const double sin_theta0 = sin(theta0);
      const double cos_theta1 = cos(theta1);
      const double sin_theta1 = sin(theta1);
      for (size_t n = 0; n < npoints - 1; n++) {
        const point_t * const p_a = points + (n + 0) % npoints;
        const point_t * const p_b = points + (n + 1) % npoints;
        const double r_a = p_a->y;
        const double z_a = p_a->x - 4.;
        const double r_b = p_b->y;
        const double z_b = p_b->x - 4.;
        const double nr_a = p_a->ny;
        const double nz_a = p_a->nx;
        const double nr_b = p_b->ny;
        const double nz_b = p_b->nx;
        const double x_a0 = r_a * cos_theta0;
        const double y_a0 = r_a * sin_theta0;
        const double x_a1 = r_a * cos_theta1;
        const double y_a1 = r_a * sin_theta1;
        const double x_b0 = r_b * cos_theta0;
        const double y_b0 = r_b * sin_theta0;
        const double x_b1 = r_b * cos_theta1;
        const double y_b1 = r_b * sin_theta1;
        const double nx_a0 = nr_a * cos_theta0;
        const double ny_a0 = nr_a * sin_theta0;
        const double nx_a1 = nr_a * cos_theta1;
        const double ny_a1 = nr_a * sin_theta1;
        const double nx_b0 = nr_b * cos_theta0;
        const double ny_b0 = nr_b * sin_theta0;
        const double nx_b1 = nr_b * cos_theta1;
        const double ny_b1 = nr_b * sin_theta1;
        render_triangle(
            &camera,
            &light,
            &screen,
            &(triangle_t){
              .vertices = {
                (vector_t){.x = x_a0, .y = y_a0, .z = z_a},
                (vector_t){.x = x_b0, .y = y_b0, .z = z_b},
                (vector_t){.x = x_a1, .y = y_a1, .z = z_a},
              },
              .vertex_normals = {
                (vector_t){.x = nx_a0, .y = ny_a0, .z = nz_a},
                (vector_t){.x = nx_b0, .y = ny_b0, .z = nz_b},
                (vector_t){.x = nx_a1, .y = ny_a1, .z = nz_a},
              },
            },
            canvas
        );
        render_triangle(
            &camera,
            &light,
            &screen,
            &(triangle_t){
              .vertices = {
                (vector_t){.x = x_b0, .y = y_b0, .z = z_b},
                (vector_t){.x = x_b1, .y = y_b1, .z = z_b},
                (vector_t){.x = x_a1, .y = y_a1, .z = z_a},
              },
              .vertex_normals = {
                (vector_t){.x = nx_b0, .y = ny_b0, .z = nz_b},
                (vector_t){.x = nx_b1, .y = ny_b1, .z = nz_b},
                (vector_t){.x = nx_a1, .y = ny_a1, .z = nz_a},
              },
            },
            canvas
        );
      }
    }
  }
  output(file_name, screen.width, screen.height, canvas);
  free(canvas);
  return 0;
}

