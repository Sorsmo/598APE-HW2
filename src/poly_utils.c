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

Poly create_poly(void) {
  Poly p;
  memset(p.coeffs, 0, sizeof(p.coeffs));
  p.degree = -1;
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
  return p.degree;
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

Poly coeff_mod(Poly p, double modulus) {
  Poly out = create_poly();
  int64_t degree = poly_degree(p);
  if (degree < 0) {
    return out;
  }
  for (int64_t i = 0; i <= degree; i++) {
    double coeff = p.coeffs[i];
    if (fabs(coeff) > COEFF_EPS) {
      double rounded = round(coeff);
      double m = positive_fmod(rounded, modulus);
      set_coeff(&out, i, m);
    }
  }
  return out;
}

Poly poly_add(Poly a, Poly b) {
  Poly sum = create_poly();
  int64_t deg_a = poly_degree(a);
  int64_t deg_b = poly_degree(b);
  int64_t max_deg = deg_a > deg_b ? deg_a : deg_b;
  if (max_deg < 0) {
    return sum;
  }
  for (int64_t i = 0; i <= max_deg; i++) {
    sum.coeffs[i] = a.coeffs[i] + b.coeffs[i];
  }
  poly_trim_to(&sum, max_deg);
  return sum;
}

Poly poly_mul_scalar(Poly p, double scalar) {
  Poly res = create_poly();
  int64_t degree = poly_degree(p);
  if (degree < 0) {
    return res;
  }
  int64_t last_nonzero = -1;
  for (int64_t i = 0; i <= degree; i++) {
    res.coeffs[i] = p.coeffs[i] * scalar;
    if (fabs(res.coeffs[i]) > COEFF_EPS) {
      last_nonzero = i;
    } else {
      res.coeffs[i] = 0.0;
    }
  }
  res.degree = last_nonzero;
  return res;
}

Poly poly_mul(Poly a, Poly b) {
  Poly res = create_poly();

  int64_t deg_a = poly_degree(a);
  int64_t deg_b = poly_degree(b);
  if (deg_a < 0 || deg_b < 0) {
    return res;
  }

  for (int64_t i = 0; i <= deg_a; i++) {
    double coeff_a = a.coeffs[i];
    if (fabs(coeff_a) > COEFF_EPS) {
      for (int64_t j = 0; j <= deg_b; j++) {
        double coeff_b = b.coeffs[j];
        if (fabs(coeff_b) > COEFF_EPS) {
          assert(i + j < MAX_POLY_DEGREE);
          res.coeffs[i + j] += coeff_a * coeff_b;
        }
      }
    }
  }
  poly_trim_to(&res, deg_a + deg_b);
  return res;
}

void poly_divmod(Poly num, Poly den, Poly *quot, Poly *rem) {
  int64_t modulus_degree = poly_degree(den);
  assert(modulus_degree > 0);
  assert(fabs(den.coeffs[0] - 1.0) <= COEFF_EPS);
  assert(fabs(den.coeffs[modulus_degree] - 1.0) <= COEFF_EPS);
  for (int64_t i = 1; i < modulus_degree; ++i) {
    assert(fabs(den.coeffs[i]) <= COEFF_EPS);
  }

  Poly remainder = num;
  if (quot != NULL) {
    *quot = create_poly();
  }

  int64_t highest_degree = poly_degree(remainder);
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

Poly poly_round_div_scalar(Poly x, double divisor) {
  Poly out = create_poly();
  assert(fabs(divisor) > 1e-9);

  int64_t degree = poly_degree(x);
  if (degree < 0) {
    return out;
  }
  for (int64_t i = 0; i <= degree; i++) {
    double v = x.coeffs[i];
    double rounded = round(v / divisor);
    if (fabs(rounded) > COEFF_EPS) {
      set_coeff(&out, i, rounded);
    }
  }
  return out;
}
