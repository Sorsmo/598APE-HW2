#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Poly create_poly(void) {
  Poly p;
  memset(p.coeffs, 0, sizeof(p.coeffs));
  return p;
}

double positive_fmod(double x, double m) {
  assert(m > 0.0);
  double r = fmod(x, m);
  if (r < 0.0)
    r += m;
  return r;
}

int64_t poly_degree(Poly p) {
  for (int64_t i = MAX_POLY_DEGREE - 1; i >= 0; i--) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      return i;
    }
  }
  return 0;
}

double get_coeff(Poly p, int64_t degree) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return 0.0;
  }
  return p.coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return;
  }
  p->coeffs[degree] = value;
}

Poly coeff_mod(Poly p, double modulus) {
  Poly out = create_poly();
  int64_t degree = poly_degree(p);
  for (int64_t i = 0; i <= degree; i++) {
    double coeff = p.coeffs[i];
    if (fabs(coeff) > 1e-9) {
      double rounded = round(coeff);
      double m = positive_fmod(rounded, modulus);
      out.coeffs[i] = m;
    }
  }
  return out;
}

Poly poly_add(Poly a, Poly b) {
  Poly sum = create_poly();
  int64_t deg_a = poly_degree(a);
  int64_t deg_b = poly_degree(b);
  int64_t max_deg = deg_a > deg_b ? deg_a : deg_b;
  for (int64_t i = 0; i <= max_deg; i++) {
    sum.coeffs[i] = a.coeffs[i] + b.coeffs[i];
  }
  return sum;
}

Poly poly_mul_scalar(Poly p, double scalar) {
  Poly res = create_poly();
  int64_t degree = poly_degree(p);
  for (int64_t i = 0; i <= degree; i++) {
    res.coeffs[i] = p.coeffs[i] * scalar;
  }
  return res;
}

Poly poly_mul(Poly a, Poly b) {
  Poly res = create_poly();

  int64_t deg_a = poly_degree(a);
  int64_t deg_b = poly_degree(b);

  for (int64_t i = 0; i <= deg_a; i++) {
    double coeff_a = a.coeffs[i];
    if (fabs(coeff_a) > 1e-9) {
      for (int64_t j = 0; j <= deg_b; j++) {
        double coeff_b = b.coeffs[j];
        if (fabs(coeff_b) > 1e-9) {
          assert(i + j < MAX_POLY_DEGREE);
          res.coeffs[i + j] += coeff_a * coeff_b;
        }
      }
    }
  }
  return res;
}

void poly_divmod(Poly num, Poly den, Poly *quot, Poly *rem) {
  const double eps = 1e-9;
  int64_t modulus_degree = poly_degree(den);
  assert(modulus_degree > 0);
  assert(fabs(den.coeffs[0] - 1.0) <= eps);
  assert(fabs(den.coeffs[modulus_degree] - 1.0) <= eps);
  for (int64_t i = 1; i < modulus_degree; ++i) {
    assert(fabs(den.coeffs[i]) <= eps);
  }

  Poly remainder = num;
  if (quot != NULL) {
    *quot = create_poly();
  }

  int64_t highest_degree = poly_degree(remainder);
  for (int64_t i = highest_degree; i >= modulus_degree; --i) {
    double coeff = remainder.coeffs[i];
    if (fabs(coeff) <= eps) {
      continue;
    }
    remainder.coeffs[i] = 0.0;
    remainder.coeffs[i - modulus_degree] -= coeff;
  }

  if (rem != NULL) {
    *rem = remainder;
  }
}

Poly poly_round_div_scalar(Poly x, double divisor) {
  Poly out = create_poly();
  assert(fabs(divisor) > 1e-9);

  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    double v = x.coeffs[i];
    out.coeffs[i] = round(v / divisor);
  }
  return out;
}
