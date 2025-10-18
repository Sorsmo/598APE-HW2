#include "he.h"
#include "poly_utils.h"
#include "ring_utils.h"
#include <math.h>
#include <stdlib.h>

double decrypt(SecretKey sk, size_t n, double q, const Poly *poly_mod, double t,
                Ciphertext ct) {
  Poly c1s;
  ring_mul_mod(&ct.c1, &sk, q, poly_mod, &c1s);
  Poly scaled_pt;
  ring_add_mod(&c1s, &ct.c0, q, poly_mod, &scaled_pt);

  Poly dec;
  poly_init(&dec);

  int64_t degree = poly_degree(&scaled_pt);
  for (int64_t i = 0; i <= degree; i++) {
    double coeff = scaled_pt.coeffs[i];
    if (fabs(coeff) > 1e-9) {
      double v = round(coeff);
      double result = round(t * v / q);
      set_coeff(&dec, i, positive_fmod(result, t));
    }
  }

  return round(get_coeff(&dec, 0));
}
