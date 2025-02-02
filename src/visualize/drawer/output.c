#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "output.h"

int output(
    const char file_name[],
    const size_t width,
    const size_t height,
    pixel_t * const canvas
) {
  const size_t nitems = width * height;
  const size_t size = sizeof(uint8_t);
  uint8_t * const buffer = malloc(nitems * size);
  if (NULL == buffer) {
    puts("failed to allocate image buffer");
    return 1;
  }
  for (size_t cnt = 0, j = 0; j < height; j++) {
    for (size_t i = 0; i < width; i++) {
      // flip in y
      const size_t index = (height - j - 1) * width + i;
      const pixel_t * const pixel = canvas + index;
      const color_t * const color = &pixel->color;
      buffer[cnt++] = color->grayscale;
    }
  }
  // dump a pgm file
  errno = 0;
  FILE * const fp = fopen(file_name, "w");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  // header
  errno = 0;
  fprintf(fp, "P5\n%zu %zu\n255\n", width, height);
  if (0 != errno) {
    perror(file_name);
    fclose(fp);
    return 1;
  }
  // contents
  const size_t retval = fwrite(buffer, size, nitems, fp);
  if (nitems != retval) {
    fprintf(stderr, "fwrite failed (%zu expected, %zu returned)\n", nitems, retval);
    fclose(fp);
    return 1;
  }
  fclose(fp);
  // clean-up
  free(buffer);
  return 0;
}

