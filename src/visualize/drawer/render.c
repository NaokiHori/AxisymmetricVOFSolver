#include <math.h>
#include "./render.h"

static inline double dblmin3 (
    const double v0,
    const double v1,
    const double v2
) {
  return fmin(v0, fmin(v1, v2));
}

static inline double dblmax3 (
    const double v0,
    const double v1,
    const double v2
) {
  return fmax(v0, fmax(v1, v2));
}

static inline int_fast32_t intmin (
    const int_fast32_t v0,
    const int_fast32_t v1
) {
  return v0 < v1 ? v0 : v1;
}

static inline int_fast32_t intmax (
    const int_fast32_t v0,
    const int_fast32_t v1
) {
  return v1 < v0 ? v0 : v1;
}

static double edge_function(
    const vector_t * const v0,
    const vector_t * const v1,
    const vector_t * const v2
) {
  return
    + (v2->x - v0->x) * (v1->y - v0->y)
    - (v2->y - v0->y) * (v1->x - v0->x);
}

// project a point p0 described on the Cartesian domain
//   to the screen
// the resulting point p1 is described on the screen coordinate system,
//   and the third element holds the depth information
static int project(
    const camera_t * const camera,
    const screen_t * const screen,
    const vector_t * const p0,
    vector_t * const p1
) {
  // vector from the camera to the given point
  // note that an arbitrary point on this line
  //   can be represented as "camera + t * ray",
  //   where "t" is a parameter
  const vector_t ray = sub(*p0, camera->position);
  // determine parameter t where ray intersects with the screen
  // a function describing the screen surface reads
  //   nx x + ny y + nz z = d
  // "t" is uniquely determined by coupling the surface and the line:
  //   normal (camera + t ray) = d
  const double nc = inner_product(screen->normal, camera->position);
  const double nr = inner_product(screen->normal, ray);
  const double d  = inner_product(screen->normal, screen->center);
  // avoid zero division
  const double small = 1.e-8;
  if (- small < nr && nr < small) {
    return 1;
  }
  const double t = (d - nc) / nr;
  // out-of-range
  if (t <= 0. || 1. <= t) {
    return 1;
  }
  // find a vector from the screen center to the intersection on the screen
  const vector_t delta = {
    + camera->position.x + t * ray.x - screen->center.x,
    + camera->position.y + t * ray.y - screen->center.y,
    + camera->position.z + t * ray.z - screen->center.z,
  };
  // the above vector components are still in the global Cartesian coordinate,
  //   which is transformed to the screen coordinate
  // this is done by computing the inner product with
  //   the horizontal and the vertical vectors respectively
  //   which are followed by appropriate normalisations
  // these vectors are NOT normalised and
  //   their l2 norms represent the lengths of the screen,
  //   which are used to normalise
  // NOTE: two screen bases are assumed to be orthogonal
  p1->x = inner_product(delta, screen->local_x) / inner_product(screen->local_x, screen->local_x);
  p1->y = inner_product(delta, screen->local_y) / inner_product(screen->local_y, screen->local_y);
  // z is used to store the depth, which is usually negative
  p1->z = nr;
  return 0;
}

int render_triangle(
    const camera_t * const camera,
    const vector_t * const light,
    const screen_t * const screen,
    const triangle_t * const triangle,
    pixel_t * const canvas
) {
  const size_t width  = screen->width;
  const size_t height = screen->height;
  const color_t facet_color = (color_t){
    .grayscale = 0xff,
  };
  // consider a triangle projected onto the screen
  // if at least one vertex returns non-zero value,
  //   I assume this triangle is out-of-range
  triangle_t projected = {0};
  for (size_t vertex_id = 0; vertex_id < 3; vertex_id++) {
    if (0 != project(
          camera,
          screen,
          triangle->vertices + vertex_id,
          projected.vertices + vertex_id
    )) {
      // this is an expected termination and thus return success
      return 0;
    }
  }
  // three vertices of the projected triangle on the screen
  // the vector components are on the screen coordinate (2D),
  //   while the z element is used to store the depth,
  //   which will be used to compare / update the z-buffer
  const vector_t * const v0 = projected.vertices + 0;
  const vector_t * const v1 = projected.vertices + 1;
  const vector_t * const v2 = projected.vertices + 2;
  // prepare bounding box
  const int_fast32_t xmin = (0.5 + dblmin3(v0->x, v1->x, v2->x)) * width  - 1;
  const int_fast32_t xmax = (0.5 + dblmax3(v0->x, v1->x, v2->x)) * width  + 1;
  const int_fast32_t ymin = (0.5 + dblmin3(v0->y, v1->y, v2->y)) * height - 1;
  const int_fast32_t ymax = (0.5 + dblmax3(v0->y, v1->y, v2->y)) * height + 1;
  // avoid out-of-bounds access
  const size_t imin = intmin( width - 1, intmax(0, xmin));
  const size_t imax = intmin( width - 1, intmax(0, xmax));
  const size_t jmin = intmin(height - 1, intmax(0, ymin));
  const size_t jmax = intmin(height - 1, intmax(0, ymax));
  // perform in-out check for each pixel inside the bounding box
  for (size_t j = jmin; j <= jmax; j++) {
    for (size_t i = imin; i <= imax; i++) {
      // target point inside the bounding box, on the screen coordinate[-0.5 : +0.5]
      const vector_t v3 = {
        1. * (i + 0.5) / width  - 0.5,
        1. * (j + 0.5) / height - 0.5,
        0., // not used
      };
      // compute weights on the barycentric coordinate
      //   using the edge function normalised by the area
      // compute inverse signed area while taking care of zero divisions
      const double area = edge_function(v0, v1, v2);
      const double area_inv = 1. / (
          area < 0. ? fmin(area, - 1.e-16)
                    : fmax(area, + 1.e-16)
      );
      const double w0 = area_inv * edge_function(v1, v2, &v3);
      const double w1 = area_inv * edge_function(v2, v0, &v3);
      const double w2 = area_inv * edge_function(v0, v1, &v3);
      // early return if this pixel is out of the triangle
      // give small negative number instead of 0 to avoid hole
      const double small = - 1.e-8;
      if (w0 < small || w1 < small || w2 < small) {
        continue;
      }
      // we know this pixel is inside the triangle now
      pixel_t * const pixel = canvas + j * width + i;
      double * const dist1 = &pixel->depth;
      // compute depth by using harmonic average in the barycentric coordinate
      const double dist0 = 1. / (
          + w0 / v0->z
          + w1 / v1->z
          + w2 / v2->z
      );
      if (dist0 < *dist1) {
        // this triangle does not come to the nearest in this pixel
        // go on to the next pixel
        continue;
      }
      // we found this facet comes the nearest
      // update the nearest distance for later elements
      *dist1 = dist0;
      // adjust facet color (make it darker) depending on
      //   the angle between the normal vector and the light
      // first obtain the local face normal
      //   by averaging three vertex normals on the barycentric coordinate
      const vector_t * const vertex_normals = triangle->vertex_normals;
      const vector_t face_normal = normalize(
          add(
            mul(w0, vertex_normals[0]),
            add(
              mul(w1, vertex_normals[1]),
              mul(w2, vertex_normals[2])
            )
          )
      );
      const double diffuse_factor = inner_product(face_normal, *light);
      const uint8_t env = 0x22;
      const double grayscale = fmax(env, diffuse_factor * facet_color.grayscale);
      pixel->color.grayscale = (uint8_t)grayscale;
    }
  }
  return 0;
}

