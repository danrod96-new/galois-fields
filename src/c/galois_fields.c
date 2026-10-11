/*
 * galois_fields.c - Arithmetic in finite (Galois) fields GF(p^n).
 *
 * See galois_fields.h for the conventions. All values are kept below
 * GF_MAX_ORDER (2^20), so products of two coefficients fit in 64 bits.
 *
 * Copyright (C) Daniel E Rodriguez - UIS.
 */
#include <stdlib.h>
#include <string.h>
#include "galois_fields.h"

/* Elements of GF(p^n) encoded as integers c0 + c1 p + ... (base p digits). */
typedef struct
{
    uint64_t p;
    int n;
    uint64_t q;
    uint32_t* alog; /* alog[i] = code of alpha^i, 0 <= i < q - 1 */
    uint32_t* log;  /* log[code] = i, for code != 0 */
} code_field;

/* ---- Integers ---------------------------------------------------------- */

int gf_is_prime(uint64_t n)
{
	uint64_t i;

	if (n < 2) {
		return 0;
	}

	if (n < 4) {
		return 1;
	}

	if (n % 2 == 0 || n % 3 == 0) {
		return 0;
	}

	for (i = 5; i <= n / i; i += 6) {
		if (n % i == 0 || n % (i + 2) == 0) {
			return 0;
		}
	}

	return 1;
}

uint64_t gf_mod(int64_t a, uint64_t m)
{
  int64_t r = a % (int64_t) m;
  
	return (uint64_t) (r < 0 ? r + (int64_t) m : r);
}

uint64_t gf_pow_mod(uint64_t b, uint64_t e, uint64_t m)
{
	uint64_t r = 1 % m;
	b %= m;

	while (e > 0) {
		if (e & 1) {
			r = (r * b) % m;
		}
		
		b = (b * b) % m;
		e >>= 1;
	}

	return r;
}

uint64_t gf_inv_mod(uint64_t a, uint64_t p)
{
  /* Fermat: a^(p-2) = a^-1 for prime p. */
  return gf_pow_mod(a % p, p - 2, p);
}

uint64_t gf_field_order(uint64_t p, int n, uint64_t limit)
{
	uint64_t q = 1;
	int i;

	for (i = 0; i < n; i++) {
		if (q > limit / p) {
			return 0;
		}

		q *= p;
	}

	return q <= limit ? q : 0;
}

/* Distinct prime factors of m, at most max of them. */
static int prime_factors(uint64_t m, uint64_t *factors, int max)
{
	int count = 0;
	uint64_t d;

	for (d = 2; d <= m / d && count < max; d++) {
		if (m % d == 0) {
			factors[count++] = d;

			while (m % d == 0) {
				m /= d;
			}
		}
	}

	if (m > 1 && count < max) {
		factors[count++] = m;
	}

	return count;
}

uint64_t gf_primitive_root(uint64_t p)
{
	uint64_t factors[64];
	int count, i, ok;
	uint64_t g;

	if (p == 2) {
			return 1;
	}

	count = prime_factors(p - 1, factors, 64);
	
	for (g = 2; g < p; g++) {
		ok = 1;
		
		for (i = 0; i < count && ok; i++) {
			if (gf_pow_mod(g, (p - 1) / factors[i], p) == 1)
			{
					ok = 0;
			}
		}

		if (ok) {
				return g;
		}
	}

	return 0;
}

/* ---- Element arithmetic ---------------------------------------------- */

/* Reduces t (len coefficients, modified in place) modulo f (degree n) and
 * writes the n remaining coefficients to r. */
static void reduce(uint64_t p, int n, const uint64_t *f, uint64_t *t, int len, uint64_t *r)
{
	uint64_t inv_lead = gf_inv_mod(f[n], p);
	uint64_t c;
	int i, j;

	for (i = len - 1; i >= n; i--) {
		c = t[i] % p;

		if (c == 0) {
			continue;
		}

		c = (c * inv_lead) % p;

		for (j = 0; j <= n; j++) {
			t[i - n + j] = (t[i - n + j] + p - (c * f[j]) % p) % p;
		}
	}

	for (j = 0; j < n; j++) {
		r[j] = j < len ? t[j] % p : 0;
	}
}

int gf_elem_is_zero(int n, const uint64_t *a)
{
	int i;

	for (i = 0; i < n; i++) {
		if (a[i] != 0) {
			return 0;
		}
	}

	return 1;
}

void gf_elem_add(uint64_t p, int n, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	int i;

	for (i = 0; i < n; i++) {
		r[i] = (a[i] + b[i]) % p;
	}
}

void gf_elem_sub(uint64_t p, int n, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	int i;

	for (i = 0; i < n; i++) {
		r[i] = (a[i] + p - b[i]) % p;
	}
}

void gf_elem_mul(uint64_t p, int n, const uint64_t *f, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t t[2 * GF_MAX_DEGREE];
	int i, j;

	memset(t, 0, sizeof(t));

	for (i = 0; i < n; i++) {
		if (a[i] == 0) {
				continue;
		}

		for (j = 0; j < n; j++) {
			t[i + j] = (t[i + j] + a[i] * b[j]) % p;
		}
	}

	reduce(p, n, f, t, 2 * n - 1, r);
}

int gf_elem_pow(uint64_t p, int n, const uint64_t *f, const uint64_t *a, uint64_t e, uint64_t *r)
{
	uint64_t base[GF_MAX_DEGREE];
	uint64_t q = gf_field_order(p, n, UINT64_MAX / 2);
	int i;

	if (gf_elem_is_zero(n, a)) {
		if (e == 0){
			return GF_ERR_DIVISION_BY_ZERO;
		}

		memset(r, 0, n * sizeof(uint64_t));
		return GF_OK;
	}

	/* Non-zero elements form a cyclic group of order q - 1. */
	if (q > 1) {
		e %= (q - 1);
	}

	memcpy(base, a, n * sizeof(uint64_t));
	memset(r, 0, n * sizeof(uint64_t));
	r[0] = 1;

	while (e > 0) {
		if (e & 1) {
			gf_elem_mul(p, n, f, r, base, r);
		}

		gf_elem_mul(p, n, f, base, base, base);
		e >>= 1;
	}

	for (i = 0; i < n; i++){
		r[i] %= p;
	}

	return GF_OK;
}

int gf_elem_div(uint64_t p, int n, const uint64_t *f, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t inverse[GF_MAX_DEGREE];
	uint64_t q = gf_field_order(p, n, UINT64_MAX / 2);

	if (gf_elem_is_zero(n, b)) {
		return GF_ERR_DIVISION_BY_ZERO;
	}

	/* b^(q - 2) = b^-1 in GF(q). */
	gf_elem_pow(p, n, f, b, q - 2, inverse);
	gf_elem_mul(p, n, f, a, inverse, r);
	
	return GF_OK;
}

/* ---- Primitive polynomials ------------------------------------------- */

static int is_one(int n, const uint64_t *a)
{
	int i;

	if (a[0] != 1) {
		return 0;
	}

	for (i = 1; i < n; i++) {
		if (a[i] != 0) {
			return 0;
		}
	}

	return 1;
}

/* f is primitive iff x has multiplicative order exactly p^n - 1 modulo f
 * (which also implies that f is irreducible). */
static int is_primitive_with_factors(const uint64_t *f, int n, uint64_t p, uint64_t q,
                                     const uint64_t *factors, int nfactors)
{
	uint64_t x_poly[GF_MAX_DEGREE + 1];
	uint64_t x[GF_MAX_DEGREE];
	uint64_t r[GF_MAX_DEGREE];
	int i;

	if (f[n] % p == 0 || f[0] % p == 0) {
		return 0;
	}

	memset(x_poly, 0, sizeof(x_poly));
	x_poly[1] = 1;
	reduce(p, n, f, x_poly, 2, x);

	/* x^(q-1) must be 1. gf_elem_pow reduces exponents modulo q - 1, so the
		* powers are computed explicitly here. */
	{
		uint64_t e = q - 1;
		uint64_t base[GF_MAX_DEGREE];

		memcpy(base, x, n * sizeof(uint64_t));
		memset(r, 0, n * sizeof(uint64_t));
		r[0] = 1;
		
		while (e > 0) {
			if (e & 1) {
				gf_elem_mul(p, n, f, r, base, r);
			}
			
			gf_elem_mul(p, n, f, base, base, base);
			e >>= 1;
		}

		if (!is_one(n, r)) {
			return 0;
		}
	}

	for (i = 0; i < nfactors; i++)
	{
		uint64_t e = (q - 1) / factors[i];
		uint64_t base[GF_MAX_DEGREE];

		memcpy(base, x, n * sizeof(uint64_t));
		memset(r, 0, n * sizeof(uint64_t));
		r[0] = 1;

		while (e > 0) {
			if (e & 1) {
				gf_elem_mul(p, n, f, r, base, r);
			}

			gf_elem_mul(p, n, f, base, base, base);
			e >>= 1;
		}

		if (is_one(n, r)) {
			return 0;
		}
	}

	return 1;
}

int gf_is_primitive_polynomial(const uint64_t *f, int n, uint64_t p)
{
	uint64_t factors[64];
	uint64_t q = gf_field_order(p, n, GF_MAX_ORDER);
	int nfactors;

	if (q == 0 || n < 1 || n > GF_MAX_DEGREE) {
		return 0;
	}

	nfactors = prime_factors(q - 1, factors, 64);

	return is_primitive_with_factors(f, n, p, q, factors, nfactors);
}

int gf_default_primitive_polynomial(uint64_t p, int n, uint64_t *f)
{
	uint64_t factors[64];
	uint64_t q, k, v;
	int nfactors, i;

	if (n == 1) {
		uint64_t g = gf_primitive_root(p);
		f[0] = (p - g) % p;
		f[1] = 1;

		return GF_OK;
	}

	q = gf_field_order(p, n, GF_MAX_ORDER);

	if (q == 0) {
		return GF_ERR_NOT_FOUND;
	}

	nfactors = prime_factors(q - 1, factors, 64);
	f[n] = 1;
	
	for (k = 1; k < q; k++) {
		v = k;

		for (i = 0; i < n; i++) {
			f[i] = v % p;
			v /= p;
		}

		if (f[0] == 0) {
			continue;
		}

		if (is_primitive_with_factors(f, n, p, q, factors, nfactors)) {
			return GF_OK;
		}
	}

	return GF_ERR_NOT_FOUND;
}

int gf_field_elements(uint64_t p, int n, const uint64_t *f, double *elements)
{
	uint64_t q = gf_field_order(p, n, GF_MAX_ORDER);
	uint64_t cur[GF_MAX_DEGREE + 1];
	uint64_t row;
	int c;

	if (q == 0) {
		return GF_ERR_NOT_FOUND;
	}

	if (n == 1)
	{
		for (row = 0; row < q; row++) {
			elements[row] = (double) row;
		}

		return GF_OK;
	}

	for (c = 0; c < n; c++){
		elements[c * q] = 0;
	}

	memset(cur, 0, sizeof(cur));
	cur[0] = 1;

	for (row = 1; row < q; row++)
	{
		for (c = 0; c < n; c++){
			elements[c * q + row] = (double)cur[c];
		}
		
		/* cur = cur * x mod f */
		for (c = n; c > 0; c--){
			cur[c] = cur[c - 1];
		}

		cur[0] = 0;
		reduce(p, n, f, cur, n + 1, cur);
		cur[n] = 0;
	}

	return GF_OK;
}

/* ---- Minimal polynomials --------------------------------------------- */

static uint64_t code_add(const code_field *F, uint64_t a, uint64_t b)
{
	uint64_t r = 0, w = 1;
	int i;

	for (i = 0; i < F->n; i++) {
		r += ((a % F->p + b % F->p) % F->p) * w;
		a /= F->p;
		b /= F->p;
		w *= F->p;
	}

	return r;
}

static uint64_t code_neg(const code_field *F, uint64_t a)
{
	uint64_t r = 0, w = 1;
	int i;

	for (i = 0; i < F->n; i++) {
		r += ((F->p - a % F->p) % F->p) * w;
		a /= F->p;
		w *= F->p;
	}

	return r;
}

static uint64_t code_mul(const code_field* F, uint64_t a, uint64_t b)
{
	if (a == 0 || b == 0) {
		return 0;
	}

	return F->alog[(F->log[a] + F->log[b]) % (F->q - 1)];
}

int gf_minimal_polynomial(uint64_t p, int n, const uint64_t *f, int64_t e, uint64_t *out)
{
	code_field F;
	uint64_t cur[GF_MAX_DEGREE + 1];
	uint64_t* m = NULL;
	uint64_t k0, k, code, w;
	int d = 0, i, c;

	memset(out, 0, (n + 1) * sizeof(uint64_t));

	if (n == 1) {
		out[0] = (p - gf_mod(e, p)) % p;
		out[1] = 1;
		return GF_OK;
	}

	if (e < 0) {
		/* Minimal polynomial of 0 is x. */
		out[1] = 1;
		return GF_OK;
	}

	F.p = p;
	F.n = n;
	F.q = gf_field_order(p, n, GF_MAX_ORDER);

	if (F.q == 0) {
		return GF_ERR_NOT_FOUND;
	}

	F.alog = (uint32_t*) malloc((F.q - 1) * sizeof(uint32_t));
	F.log = (uint32_t*) calloc(F.q, sizeof(uint32_t));
	m = (uint64_t*) calloc(n + 1, sizeof(uint64_t));

	if (F.alog == NULL || F.log == NULL || m == NULL) {
		free(F.alog);
		free(F.log);
		free(m);
		return GF_ERR_MEMORY;
	}

	/* Power table alpha^i, alpha a root of f. */
	memset(cur, 0, sizeof(cur));
	cur[0] = 1;

	for (k = 0; k < F.q - 1; k++) {
		code = 0;
		w = 1;

		for (c = 0; c < n; c++) {
			code += cur[c] * w;
			w *= p;
		}

		F.alog[k] = (uint32_t)code;
		F.log[code] = (uint32_t)k;

		for (c = n; c > 0; c--){
			cur[c] = cur[c - 1];
		}

		cur[0] = 0;
		reduce(p, n, f, cur, n + 1, cur);
		cur[n] = 0;
	}

	/* m(x) = product of (x - alpha^(e p^j)) over the distinct conjugates. */
	m[0] = 1;
	k0 = (uint64_t)e % (F.q - 1);
	k = k0;

	do {
		uint64_t root = F.alog[k];
		uint64_t minus_root = code_neg(&F, root);

		/* m = m * (x - root) */
		m[d + 1] = m[d];

		for (i = d; i > 0; i--) {
			m[i] = code_add(&F, m[i - 1], code_mul(&F, minus_root, m[i]));
		}

		m[0] = code_mul(&F, minus_root, m[0]);
		d++;
		k = (k * p) % (F.q - 1);
	} while (k != k0 && d < n);

	/* The coefficients lie in the prime subfield, so their codes are < p. */
	for (i = 0; i <= d; i++) {
		out[i] = m[i] % p;
	}

	free(F.alog);
	free(F.log);
	free(m);

	return GF_OK;
}

/* ---- Polynomials over GF(p) ------------------------------------------ */

int gf_poly_degree(const uint64_t *a, int len)
{
	int i;

	for (i = len - 1; i >= 0; i--){
		if (a[i] != 0) {
			return i;
		}
	}

	return -1;
}

void gf_poly_add(uint64_t p, const uint64_t *a, int la, const uint64_t *b, int lb, uint64_t *r)
{
	int i, len = la > lb ? la : lb;
	
	for (i = 0; i < len; i++) {
		r[i] = ((i < la ? a[i] : 0) + (i < lb ? b[i] : 0)) % p;
	}
}

void gf_poly_sub(uint64_t p, const uint64_t *a, int la, const uint64_t *b, int lb, uint64_t *r)
{
	int i, len = la > lb ? la : lb;

	for (i = 0; i < len; i++) {
		r[i] = ((i < la ? a[i] : 0) + p - (i < lb ? b[i] : 0)) % p;
	}
}

void gf_poly_mul(uint64_t p, const uint64_t *a, int la, const uint64_t *b, int lb, uint64_t *r)
{
	int i, j;
	memset(r, 0, (la + lb - 1) * sizeof(uint64_t));

	for (i = 0; i < la; i++) {
		if (a[i] == 0) {
			continue;
		}

		for (j = 0; j < lb; j++) {
			r[i + j] = (r[i + j] + a[i] * b[j]) % p;
		}
	}
}

int gf_poly_divmod(uint64_t p, const uint64_t *a, int la, const uint64_t *b, int lb, uint64_t *q, uint64_t *r)
{
	int db = gf_poly_degree(b, lb);
	int i, j;
	uint64_t inv_lead, c;

	if (db < 0) {
		return GF_ERR_DIVISION_BY_ZERO;
	}

	memcpy(r, a, la * sizeof(uint64_t));
	memset(q, 0, la * sizeof(uint64_t));
	inv_lead = gf_inv_mod(b[db], p);
	
	for (i = la - 1; i >= db; i--) {
		c = (r[i] * inv_lead) % p;
		
		if (c == 0) {
			continue;
		}

		q[i - db] = c;
		
		for (j = 0; j <= db; j++) {
			r[i - db + j] = (r[i - db + j] + p - (c * b[j]) % p) % p;
		}
	}

	return GF_OK;
}
