#if !defined(VECTOR_H)
#define VECTOR_H

typedef struct {
  double x;
  double y;
  double z;
} vector_t;

extern vector_t normalize(
    const vector_t v
);

extern vector_t add(
    const vector_t v0,
    const vector_t v1
);

extern vector_t sub(
    const vector_t v0,
    const vector_t v1
);

extern vector_t mul(
    const double c,
    const vector_t v
);

extern double inner_product(
    const vector_t v0,
    const vector_t v1
);

extern vector_t outer_product(
    const vector_t v0,
    const vector_t v1
);

#endif // VECTOR_H
