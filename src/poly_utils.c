#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

static const double COEFF_EPS = 1e-9;

static void poly_trim_to(Poly *p, int64_t max_index) {
  if (max_index < 0) {
    p->degree = -1;
    return;
  }
  if (max_index >= MAX_POLY_DEGREE) {
    max_index = MAX_POLY_DEGREE - 1;
  }
  for (int64_t i = max_index; i >= 0; --i) {
    if (fabs(p->coeffs[i]) > COEFF_EPS) {
      p->degree = i;
      return;
    }
    p->coeffs[i] = 0.0;
  }
  p->degree = -1;
}

void poly_init(Poly *p) {
  if (p == NULL) {
    return;
  }
  memset(p->coeffs, 0, sizeof(p->coeffs));
  p->degree = -1;
}

double positive_fmod(double x, double m) {
  assert(m > 0.0);
  double r = fmod(x, m);
  if (r < 0.0)
    r += m;
  return r;
}

int64_t poly_degree(const Poly *p) {
  if (p == NULL) {
    return -1;
  }
  return p->degree;
}

double get_coeff(const Poly *p, int64_t degree) {
  if (p == NULL || degree >= MAX_POLY_DEGREE || degree < 0) {
    return 0.0;
  }
  return p->coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (p == NULL || degree >= MAX_POLY_DEGREE || degree < 0) {
    return;
  }
  if (fabs(value) <= COEFF_EPS) {
    p->coeffs[degree] = 0.0;
    if (p->degree == degree) {
      poly_trim_to(p, degree - 1);
    }
    return;
  }

  p->coeffs[degree] = value;
  if (degree > p->degree) {
    p->degree = degree;
  }
}

void coeff_mod(const Poly *p, double modulus, Poly *out) {
  if (out == NULL) {
    return;
  }
  poly_init(out);
  int64_t degree = poly_degree(p);
  if (degree < 0) {
    return;
  }
  for (int64_t i = 0; i <= degree; i++) {
    double coeff = p->coeffs[i];
    if (fabs(coeff) > COEFF_EPS) {
      double rounded = round(coeff);
      double m = positive_fmod(rounded, modulus);
      set_coeff(out, i, m);
    }
  }
}

void poly_add(const Poly *a, const Poly *b, Poly *out) {
  if (out == NULL) {
    return;
  }
  poly_init(out);
  int64_t deg_a = poly_degree(a);
  int64_t deg_b = poly_degree(b);
  int64_t max_deg = deg_a > deg_b ? deg_a : deg_b;
  if (max_deg < 0) {
    return;
  }
  for (int64_t i = 0; i <= max_deg; i++) {
    double av = (a && i <= deg_a) ? a->coeffs[i] : 0.0;
    double bv = (b && i <= deg_b) ? b->coeffs[i] : 0.0;
    out->coeffs[i] = av + bv;
  }
  poly_trim_to(out, max_deg);
}

void poly_mul_scalar(const Poly *p, double scalar, Poly *out) {
  if (out == NULL) {
    return;
  }
  poly_init(out);
  int64_t degree = poly_degree(p);
  if (degree < 0) {
    return;
  }
  int64_t last_nonzero = -1;
  for (int64_t i = 0; i <= degree; i++) {
    double coeff = p->coeffs[i] * scalar;
    if (fabs(coeff) > COEFF_EPS) {
      out->coeffs[i] = coeff;
      last_nonzero = i;
    }
  }
  out->degree = last_nonzero;
}

void poly_mul(const Poly *a, const Poly *b, Poly *out) {
  if (out == NULL) {
    return;
  }
  poly_init(out);

  int64_t deg_a = poly_degree(a);
  int64_t deg_b = poly_degree(b);
  if (deg_a < 0 || deg_b < 0) {
    return;
  }

  for (int64_t i = 0; i <= deg_a; i++) {
    double coeff_a = a->coeffs[i];
    if (fabs(coeff_a) > COEFF_EPS) {
      for (int64_t j = 0; j <= deg_b; j++) {
        double coeff_b = b->coeffs[j];
        if (fabs(coeff_b) > COEFF_EPS) {
          assert(i + j < MAX_POLY_DEGREE);
          out->coeffs[i + j] += coeff_a * coeff_b;
        }
      }
    }
  }
  poly_trim_to(out, deg_a + deg_b);
}

void poly_divmod(const Poly *numerator, const Poly *denominator, Poly *quot,
                 Poly *rem) {
  int64_t modulus_degree = poly_degree(denominator);
  assert(modulus_degree > 0);
  assert(fabs(denominator->coeffs[0] - 1.0) <= COEFF_EPS);
  assert(fabs(denominator->coeffs[modulus_degree] - 1.0) <= COEFF_EPS);
  for (int64_t i = 1; i < modulus_degree; ++i) {
    assert(fabs(denominator->coeffs[i]) <= COEFF_EPS);
  }

  Poly remainder = *numerator;
  if (quot != NULL) {
    poly_init(quot);
  }

  int64_t highest_degree = poly_degree(&remainder);
  if (highest_degree < 0) {
    if (rem != NULL) {
      *rem = remainder;
    }
    return;
  }
  for (int64_t i = highest_degree; i >= modulus_degree; --i) {
    double coeff = remainder.coeffs[i];
    if (fabs(coeff) <= COEFF_EPS) {
      continue;
    }
    remainder.coeffs[i] = 0.0;
    remainder.coeffs[i - modulus_degree] -= coeff;
  }
  poly_trim_to(&remainder, modulus_degree - 1);

  if (rem != NULL) {
    *rem = remainder;
  }
}

void poly_round_div_scalar(const Poly *x, double divisor, Poly *out) {
  if (out == NULL) {
    return;
  }
  poly_init(out);
  assert(fabs(divisor) > 1e-9);

  int64_t degree = poly_degree(x);
  if (degree < 0) {
    return;
  }
  for (int64_t i = 0; i <= degree; i++) {
    double v = x->coeffs[i];
    double rounded = round(v / divisor);
    if (fabs(rounded) > COEFF_EPS) {
      set_coeff(out, i, rounded);
    }
  }
}
