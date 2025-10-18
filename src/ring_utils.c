#include "ring_utils.h"
#include "poly_utils.h"

void ring_add_mod(const Poly *x, const Poly *y, double modulus,
                  const Poly *poly_mod, Poly *out) {
  Poly sum;
  poly_add(x, y, &sum);

  Poly sum_mod;
  coeff_mod(&sum, modulus, &sum_mod);

  Poly rem;
  poly_divmod(&sum_mod, poly_mod, NULL, &rem);

  coeff_mod(&rem, modulus, out);
}

void ring_mul_mod(const Poly *x, const Poly *y, double modulus,
                  const Poly *poly_mod, Poly *out) {
  Poly prod;
  poly_mul(x, y, &prod);

  Poly prod_mod;
  coeff_mod(&prod, modulus, &prod_mod);

  Poly rem;
  poly_divmod(&prod_mod, poly_mod, NULL, &rem);

  coeff_mod(&rem, modulus, out);
}

void ring_mul_no_mod_q(const Poly *x, const Poly *y, const Poly *poly_mod,
                       Poly *out) {
  Poly prod;
  poly_mul(x, y, &prod);
  poly_divmod(&prod, poly_mod, NULL, out);
}

void ring_add_no_mod_q(const Poly *x, const Poly *y, const Poly *poly_mod,
                       Poly *out) {
  Poly sum;
  poly_add(x, y, &sum);
  poly_divmod(&sum, poly_mod, NULL, out);
}

void ring_mul_poly_mod(const Poly *x, const Poly *y, const Poly *poly_mod,
                       Poly *out) {
  Poly prod;
  poly_mul(x, y, &prod);
  poly_divmod(&prod, poly_mod, NULL, out);
}

void ring_add_poly_mod(const Poly *x, const Poly *y, const Poly *poly_mod,
                       Poly *out) {
  Poly sum;
  poly_add(x, y, &sum);
  poly_divmod(&sum, poly_mod, NULL, out);
}
