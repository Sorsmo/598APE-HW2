#ifndef HE_H
#define HE_H

#include "types.h"
#include <stdint.h>

typedef struct {
  PublicKey pk;
  SecretKey sk;
} KeyPair;

KeyPair keygen(size_t n, double q, const Poly *poly_mod);

Ciphertext encrypt(PublicKey pk, size_t n, double q, const Poly *poly_mod, double t,
                   double pt);

double decrypt(SecretKey sk, size_t n, double q, const Poly *poly_mod, double t,
               Ciphertext ct);

void encode_plain_integer(double t, double pt, Poly *out);

Ciphertext add_plain(Ciphertext ct, double q, double t, const Poly *poly_mod,
                     double pt);

Ciphertext add_cipher(Ciphertext c1, Ciphertext c2, double q, const Poly *poly_mod);

Ciphertext mul_plain(Ciphertext ct, double q, double t, const Poly *poly_mod,
                     double pt);

EvalKey evaluate_keygen(SecretKey sk, size_t n, double q, const Poly *poly_mod,
                        double p);

Ciphertext mul_cipher(Ciphertext c1, Ciphertext c2, double q, double t,
                      double p, const Poly *poly_mod, EvalKey rlk);

#endif
