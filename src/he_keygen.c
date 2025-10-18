#include "he.h"
#include "poly_random.h"
#include "poly_utils.h"
#include "ring_utils.h"

KeyPair keygen(size_t n, double q, const Poly *poly_mod) {
  SecretKey s = gen_binary_poly(n);
  Poly a = gen_uniform_poly(n, q);
  Poly e = gen_normal_poly(n, 0.0, 1.0);

  Poly neg_a;
  poly_mul_scalar(&a, -1, &neg_a);
  Poly as;
  ring_mul_mod(&neg_a, &s, q, poly_mod, &as);
  Poly neg_e;
  poly_mul_scalar(&e, -1, &neg_e);
  Poly b;
  ring_add_mod(&as, &neg_e, q, poly_mod, &b);

  KeyPair keys;
  keys.pk.a = a;
  keys.pk.b = b;
  keys.sk = s;

  return keys;
}

EvalKey evaluate_keygen(SecretKey sk, size_t n, double q, const Poly *poly_mod,
                        double p) {
  double new_modulus = q * p;
  Poly a = gen_uniform_poly(n, new_modulus);
  Poly e = gen_normal_poly(n, 0.0, 1.0);

  Poly s2;
  poly_mul(&sk, &sk, &s2);
  Poly secret_scaled;
  poly_mul_scalar(&s2, p, &secret_scaled);

  Poly neg_a;
  poly_mul_scalar(&a, -1.0, &neg_a);
  Poly as;
  ring_mul_no_mod_q(&neg_a, &sk, poly_mod, &as);
  Poly neg_e;
  poly_mul_scalar(&e, -1.0, &neg_e);
  Poly as_nege;
  ring_add_no_mod_q(&as, &neg_e, poly_mod, &as_nege);
  Poly b_ring;
  ring_add_no_mod_q(&as_nege, &secret_scaled, poly_mod, &b_ring);

  Poly b;
  coeff_mod(&b_ring, new_modulus, &b);

  EvalKey rlk;
  rlk.a = a;
  rlk.b = b;
  return rlk;
}
