/*
 * galois_fields.h - Arithmetic in finite (Galois) fields GF(p^n).
 *
 * Conventions used by the whole library:
 *  - Polynomials and extension field elements are arrays of coefficients in
 *    ascending order: a[0] is the constant term, a[len - 1] the leading one.
 *  - Elements of GF(p) are integers 0 .. p - 1.
 *  - Elements of GF(p^n), n > 1, are polynomials of degree < n reduced
 *    modulo the field's primitive polynomial, stored as n coefficients.
 *  - The element table of GF(p^n), n > 1, lists 0 first and then
 *    alpha^0, alpha^1, ..., alpha^(p^n - 2), alpha being a root of the
 *    primitive polynomial (exponential order). For n = 1 it lists 0 .. p - 1.
 *
 * Copyright (C) Daniel E Rodriguez - UIS.
 */
#ifndef GALOIS_FIELDS_H
#define GALOIS_FIELDS_H

#include <stdint.h>

#ifdef _MSC_VER
#define GF_EXPORT __declspec(dllexport)
#else
#define GF_EXPORT
#endif

/* Largest field order p^n handled (element tables have p^n rows). */
#define GF_MAX_ORDER 1048576ULL
/* Largest extension degree n. */
#define GF_MAX_DEGREE 20

/* Status codes. */
#define GF_OK 0
#define GF_ERR_MEMORY 1
#define GF_ERR_DIVISION_BY_ZERO 2
#define GF_ERR_NOT_FOUND 3

/* ---- Integers ---------------------------------------------------------- */

/* 1 if n is prime, 0 otherwise. */
GF_EXPORT int gf_is_prime(uint64_t n);

/* Non-negative remainder of a modulo m (m > 0). */
GF_EXPORT uint64_t gf_mod(int64_t a, uint64_t m);

/* b^e mod m. */
GF_EXPORT uint64_t gf_pow_mod(uint64_t b, uint64_t e, uint64_t m);

/* Inverse of a (non-zero) modulo the prime p. */
GF_EXPORT uint64_t gf_inv_mod(uint64_t a, uint64_t p);

/* p^n if it does not exceed limit, 0 otherwise. */
GF_EXPORT uint64_t gf_field_order(uint64_t p, int n, uint64_t limit);

/* Smallest primitive root modulo the prime p. */
GF_EXPORT uint64_t gf_primitive_root(uint64_t p);

/* ---- Primitive polynomials and fields ---------------------------------- */

/* 1 if f (n + 1 coefficients, f[n] != 0) is a primitive polynomial of
 * degree n over GF(p). p^n must not exceed GF_MAX_ORDER. */
GF_EXPORT int gf_is_primitive_polynomial(const uint64_t* f, int n, uint64_t p);

/* Writes the default primitive polynomial of degree n over GF(p) to f
 * (n + 1 coefficients): x - g for n = 1, g the smallest primitive root;
 * otherwise the first monic primitive polynomial when the coefficients
 * f[0] + f[1] p + ... + f[n-1] p^(n-1) are enumerated in increasing order. */
GF_EXPORT int gf_default_primitive_polynomial(uint64_t p, int n, uint64_t* f);

/* Fills the element table of GF(p^n). elements has p^n * n doubles, stored
 * column-major (Scilab layout): row r, coefficient c at elements[c * q + r]. */
GF_EXPORT int gf_field_elements(uint64_t p, int n, const uint64_t* f, double* elements);

/* ---- Element arithmetic (a, b, r: n coefficients each) ---------------- */

GF_EXPORT void gf_elem_add(uint64_t p, int n, const uint64_t* a, const uint64_t* b, uint64_t* r);
GF_EXPORT void gf_elem_sub(uint64_t p, int n, const uint64_t* a, const uint64_t* b, uint64_t* r);
GF_EXPORT void gf_elem_mul(uint64_t p, int n, const uint64_t* f, const uint64_t* a, const uint64_t* b, uint64_t* r);

/* GF_ERR_DIVISION_BY_ZERO if b is zero. */
GF_EXPORT int gf_elem_div(uint64_t p, int n, const uint64_t* f, const uint64_t* a, const uint64_t* b, uint64_t* r);

/* a^e. GF_ERR_DIVISION_BY_ZERO for 0^0, which is undefined here. */
GF_EXPORT int gf_elem_pow(uint64_t p, int n, const uint64_t* f, const uint64_t* a, uint64_t e, uint64_t* r);
GF_EXPORT int gf_elem_is_zero(int n, const uint64_t* a);

/* Minimal polynomial over GF(p) of alpha^e (n > 1, e = -1 meaning 0), or of
 * the integer e (n = 1). out receives n + 1 coefficients, zero padded. */
GF_EXPORT int gf_minimal_polynomial(uint64_t p, int n, const uint64_t* f, int64_t e, uint64_t* out);

/* ---- Polynomials over GF(p) (any length, ascending order) ------------- */

/* Degree of a (len coefficients), -1 for the zero polynomial. */
GF_EXPORT int gf_poly_degree(const uint64_t* a, int len);

/* r = a + b and r = a - b; r needs max(la, lb) coefficients. */
GF_EXPORT void gf_poly_add(uint64_t p, const uint64_t* a, int la, const uint64_t* b, int lb, uint64_t* r);
GF_EXPORT void gf_poly_sub(uint64_t p, const uint64_t* a, int la, const uint64_t* b, int lb, uint64_t* r);

/* r = a * b; r needs la + lb - 1 coefficients. */
GF_EXPORT void gf_poly_mul(uint64_t p, const uint64_t* a, int la, const uint64_t* b, int lb, uint64_t* r);

/* a = q * b + r with deg r < deg b. q needs la coefficients, r needs la.
 * GF_ERR_DIVISION_BY_ZERO if b is zero. */
GF_EXPORT int gf_poly_divmod(uint64_t p, const uint64_t* a, int la, const uint64_t* b, int lb, uint64_t* q, uint64_t* r);

#endif /* GALOIS_FIELDS_H */
