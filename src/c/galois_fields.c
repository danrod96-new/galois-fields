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
struct gf_tables {
	uint64_t p;
	int n;
	uint64_t q;
	uint64_t w[GF_MAX_DEGREE]; /* w[i] = p^i */
	uint32_t *alog; /* alog[i] = code of alpha^i, 0 <= i < q - 1 */
	uint32_t *log;  /* log[code] = i, for code != 0 */
};

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
	/* Extended Euclid: keeps t * a = r (mod p) for both rows. */
	int64_t t0 = 0, t1 = 1, t;
	uint64_t r0 = p, r1 = a % p, k, r;

	while (r1 != 0) {
		k = r0 / r1;
		r = r0 - k * r1;
		r0 = r1;
		r1 = r;
		t = t0 - (int64_t) k * t1;
		t0 = t1;
		t1 = t;
	}

	return (uint64_t) (t0 < 0 ? t0 + (int64_t) p : t0);
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
			if (gf_pow_mod(g, (p - 1) / factors[i], p) == 1) {
				ok = 0;
			}
		}

		if (ok) {
			return g;
		}
	}

	return 0;
}

/* Reduces t (len coefficients, modified in place) modulo f (degree n) and
 * writes the n remaining coefficients to r. */
static void reduce(uint64_t p, int n, const uint64_t *f, uint64_t *t, int len, uint64_t *r)
{
	/* Primitive polynomials are usually monic, which needs no inverse. */
	uint64_t inv_lead = f[n] == 1 ? 1 : gf_inv_mod(f[n], p);
	uint64_t c;
	int i, j;

	for (i = len - 1; i >= n; i--) {
		c = t[i] % p;

		if (c == 0) {
			continue;
		}

		if (inv_lead != 1) {
			c = (c * inv_lead) % p;
		}

		/* t -= c x^(i-n) f, written as t += (p - c) x^(i-n) f. */
		c = p - c;

		for (j = 0; j < n; j++) {
			t[i - n + j] = (t[i - n + j] + c * f[j]) % p;
		}

		t[i] = 0;
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

/* GF(2^n): an element packed into a word, bit i holding coefficient i. */
static uint32_t gf2_pack(int n, const uint64_t *a)
{
	uint32_t x = 0;
	int i;

	for (i = 0; i < n; i++) {
		x |= (uint32_t) (a[i] & 1) << i;
	}

	return x;
}

static void gf2_unpack(int n, uint32_t x, uint64_t *r)
{
	int i;

	for (i = 0; i < n; i++) {
		r[i] = (x >> i) & 1;
	}
}

/* a * b modulo poly (the packed primitive polynomial, degree n). */
static uint32_t gf2_mul(int n, uint32_t poly, uint32_t a, uint32_t b)
{
	uint32_t r = 0, top = (uint32_t) 1 << n;

	while (b != 0) {
		if (b & 1) {
			r ^= a;
		}

		b >>= 1;
		a <<= 1;

		if (a & top) {
			a ^= poly;
		}
	}

	return r;
}

static int gf2_degree(uint32_t x)
{
	int d = -1;

	while (x != 0) {
		x >>= 1;
		d++;
	}

	return d;
}

/* Inverse of a non-zero a modulo the irreducible poly (extended Euclid). */
static uint32_t gf2_inv(uint32_t poly, uint32_t a)
{
	uint32_t u = a, v = poly, g1 = 1, g2 = 0, tmp;
	int j;

	while (u != 1) {
		j = gf2_degree(u) - gf2_degree(v);

		if (j < 0) {
			tmp = u; u = v; v = tmp;
			tmp = g1; g1 = g2; g2 = tmp;
			j = -j;
		}

		u ^= v << j;
		g1 ^= g2 << j;
	}

	return g1;
}

static uint32_t gf2_pow(int n, uint32_t poly, uint32_t a, uint64_t e)
{
	uint32_t r = 1;

	while (e > 0) {
		if (e & 1) {
			r = gf2_mul(n, poly, r, a);
		}

		a = gf2_mul(n, poly, a, a);
		e >>= 1;
	}

	return r;
}

/* Coefficients are in 0 .. p - 1, so sums need one conditional subtraction. */
void gf_elem_add(uint64_t p, int n, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t s;
	int i;

	for (i = 0; i < n; i++) {
		s = a[i] + b[i];
		r[i] = s >= p ? s - p : s;
	}
}

void gf_elem_sub(uint64_t p, int n, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	int i;

	for (i = 0; i < n; i++) {
		r[i] = a[i] >= b[i] ? a[i] - b[i] : a[i] + p - b[i];
	}
}

void gf_elem_mul(uint64_t p, int n, const uint64_t *f, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t t[2 * GF_MAX_DEGREE];
	int i, j;

	if (n == 1) {
		r[0] = (a[0] * b[0]) % p;
		return;
	}

	if (p == 2) {
		gf2_unpack(n, gf2_mul(n, gf2_pack(n + 1, f), gf2_pack(n, a), gf2_pack(n, b)), r);
		return;
	}

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

	if (n == 1) {
		r[0] = gf_pow_mod(a[0], e, p);
		return GF_OK;
	}

	if (p == 2) {
		gf2_unpack(n, gf2_pow(n, gf2_pack(n + 1, f), gf2_pack(n, a), e), r);
		return GF_OK;
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

/* Inverse of a non-zero b modulo f (degree n) over GF(p), by the extended
 * Euclidean algorithm on polynomials: s * b = r (mod f) holds for both rows.
 * Returns 0 if gcd(b, f) is not a constant (f not irreducible). */
static int poly_inv(uint64_t p, int n, const uint64_t *f, const uint64_t *b, uint64_t *inv)
{
	uint64_t ra[GF_MAX_DEGREE + 1], rb[GF_MAX_DEGREE + 1];
	uint64_t sa[2 * GF_MAX_DEGREE + 2], sb[2 * GF_MAX_DEGREE + 2];
	uint64_t *r0 = ra, *r1 = rb, *s0 = sa, *s1 = sb, *tp, c;
	int d0 = n, d1, l0 = 0, l1 = 1, shift, j, ti;

	memcpy(r0, f, (n + 1) * sizeof(uint64_t));
	memcpy(r1, b, n * sizeof(uint64_t));
	r1[n] = 0;
	memset(sa, 0, sizeof(sa));
	memset(sb, 0, sizeof(sb));
	s1[0] = 1;
	d1 = gf_poly_degree(r1, n);

	while (d0 > 0 && d1 > 0) {
		if (d0 < d1) {
			tp = r0; r0 = r1; r1 = tp;
			tp = s0; s0 = s1; s1 = tp;
			ti = d0; d0 = d1; d1 = ti;
			ti = l0; l0 = l1; l1 = ti;
		}

		/* r0 -= c x^shift r1 and s0 -= c x^shift s1, cancelling r0's lead. */
		shift = d0 - d1;
		c = p - (r0[d0] * gf_inv_mod(r1[d1], p)) % p;

		for (j = 0; j <= d1; j++) {
			r0[shift + j] = (r0[shift + j] + c * r1[j]) % p;
		}

		if (l1 + shift > 2 * GF_MAX_DEGREE + 2) {
			return 0;
		}

		for (j = 0; j < l1; j++) {
			s0[shift + j] = (s0[shift + j] + c * s1[j]) % p;
		}

		if (l1 + shift > l0) {
			l0 = l1 + shift;
		}

		d0 = gf_poly_degree(r0, d0);
	}

	if (d0 != 0) {
		if (d1 != 0) {
			return 0;
		}

		r0 = r1;
		s0 = s1;
		l0 = l1;
	}

	/* r0 is a non-zero constant: inv = s0 / r0[0]. */
	c = gf_inv_mod(r0[0], p);

	for (j = 0; j < l0; j++) {
		s0[j] = (s0[j] * c) % p;
	}

	reduce(p, n, f, s0, l0, inv);

	return 1;
}

int gf_elem_div(uint64_t p, int n, const uint64_t *f, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t inverse[GF_MAX_DEGREE];
	uint64_t q = gf_field_order(p, n, UINT64_MAX / 2);
	uint32_t poly;

	if (gf_elem_is_zero(n, b)) {
		return GF_ERR_DIVISION_BY_ZERO;
	}

	if (n == 1) {
		r[0] = (a[0] * gf_inv_mod(b[0], p)) % p;
		return GF_OK;
	}

	if (p == 2) {
		poly = gf2_pack(n + 1, f);
		gf2_unpack(n, gf2_mul(n, poly, gf2_pack(n, a), gf2_inv(poly, gf2_pack(n, b))), r);
		return GF_OK;
	}

	/* gf_tab_div is faster still when tables are available. If f is not
	 * irreducible, fall back to b^(q - 2) as before. */
	if (!poly_inv(p, n, f, b, inverse)) {
		gf_elem_pow(p, n, f, b, q - 2, inverse);
	}

	gf_elem_mul(p, n, f, a, inverse, r);
	
	return GF_OK;
}

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

/* f is primitive if x has multiplicative order exactly p^n - 1 modulo f
 * (which also implies that f is irreducible). */
static int is_primitive_with_factors(const uint64_t *f, int n, uint64_t p, uint64_t q, const uint64_t *factors, int nfactors)
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

	if (n == 1) {
		for (row = 0; row < q; row++) {
			elements[row] = (double) row;
		}

		return GF_OK;
	}

	for (c = 0; c < n; c++){
		elements[c * q] = 0;
	}

	if (p == 2) {
		/* Bit-packed: multiplying by x is a shift, reducing it an XOR. */
		uint32_t poly = gf2_pack(n + 1, f), x = 1, top = (uint32_t) 1 << n;

		for (row = 1; row < q; row++) {
			for (c = 0; c < n; c++){
				elements[c * q + row] = (double) ((x >> c) & 1);
			}

			x <<= 1;

			if (x & top) {
				x ^= poly;
			}
		}

		return GF_OK;
	}

	memset(cur, 0, sizeof(cur));
	cur[0] = 1;

	for (row = 1; row < q; row++) {
		for (c = 0; c < n; c++){
			elements[c * q + row] = (double) cur[c];
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

gf_tables *gf_tables_create(uint64_t p, int n, const uint64_t *f)
{
	gf_tables *T;
	uint64_t cur[GF_MAX_DEGREE + 1];
	uint64_t k, code, w;
	int c;

	T = (gf_tables *) calloc(1, sizeof(gf_tables));

	if (T == NULL) {
		return NULL;
	}

	T->p = p;
	T->n = n;
	T->q = gf_field_order(p, n, GF_MAX_ORDER);

	if (T->q < 2) {
		free(T);
		return NULL;
	}

	T->alog = (uint32_t *) malloc((T->q - 1) * sizeof(uint32_t));
	T->log = (uint32_t *) calloc(T->q, sizeof(uint32_t));

	if (T->alog == NULL || T->log == NULL) {
		gf_tables_free(T);
		return NULL;
	}

	for (c = 0, w = 1; c < n; c++, w *= p) {
		T->w[c] = w;
	}

	if (n == 1) {
		/* Powers of the root of f = x - g, which is g. */
		uint64_t g = (p - f[0] % p) * gf_inv_mod(f[1], p) % p;

		for (k = 0, code = 1; k < T->q - 1; k++) {
			T->alog[k] = (uint32_t) code;
			T->log[code] = (uint32_t) k;
			code = code * g % p;
		}

		return T;
	}

	/* Power table alpha^i, alpha a root of f: cur = cur * x mod f. */
	memset(cur, 0, sizeof(cur));
	cur[0] = 1;

	for (k = 0; k < T->q - 1; k++) {
		code = 0;

		for (c = 0; c < n; c++) {
			code += cur[c] * T->w[c];
		}

		T->alog[k] = (uint32_t) code;
		T->log[code] = (uint32_t) k;

		for (c = n; c > 0; c--){
			cur[c] = cur[c - 1];
		}

		cur[0] = 0;
		reduce(p, n, f, cur, n + 1, cur);
		cur[n] = 0;
	}

	return T;
}

void gf_tables_free(gf_tables *T)
{
	if (T != NULL) {
		free(T->alog);
		free(T->log);
		free(T);
	}
}

static uint64_t tab_encode(const gf_tables *T, const uint64_t *a)
{
	uint64_t code = 0;
	int i;

	for (i = 0; i < T->n; i++) {
		code += a[i] * T->w[i];
	}

	return code;
}

static void tab_decode(const gf_tables *T, uint64_t code, uint64_t *r)
{
	int i;

	for (i = 0; i < T->n; i++) {
		r[i] = code % T->p;
		code /= T->p;
	}
}

void gf_tab_mul(const gf_tables *T, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t x = tab_encode(T, a), y = tab_encode(T, b), s;

	if (x == 0 || y == 0) {
		memset(r, 0, T->n * sizeof(uint64_t));
		return;
	}

	s = T->log[x] + T->log[y];
	tab_decode(T, T->alog[s >= T->q - 1 ? s - (T->q - 1) : s], r);
}

int gf_tab_div(const gf_tables *T, const uint64_t *a, const uint64_t *b, uint64_t *r)
{
	uint64_t x = tab_encode(T, a), y = tab_encode(T, b), s;

	if (y == 0) {
		return GF_ERR_DIVISION_BY_ZERO;
	}

	if (x == 0) {
		memset(r, 0, T->n * sizeof(uint64_t));
		return GF_OK;
	}

	s = T->log[x] + (T->q - 1) - T->log[y];
	tab_decode(T, T->alog[s >= T->q - 1 ? s - (T->q - 1) : s], r);

	return GF_OK;
}

int gf_tab_pow(const gf_tables *T, const uint64_t *a, uint64_t e, uint64_t *r)
{
	uint64_t x = tab_encode(T, a);

	if (x == 0) {
		if (e == 0) {
			return GF_ERR_DIVISION_BY_ZERO;
		}

		memset(r, 0, T->n * sizeof(uint64_t));
		return GF_OK;
	}

	tab_decode(T, T->alog[(T->log[x] * (e % (T->q - 1))) % (T->q - 1)], r);

	return GF_OK;
}

int gf_minimal_polynomial(uint64_t p, int n, const uint64_t *f, int64_t e, uint64_t *out)
{
	/* m(x) has coefficients in GF(p^n): coefficient i is m[i * n .. i * n + n - 1]. */
	uint64_t m[(GF_MAX_DEGREE + 1) * GF_MAX_DEGREE];
	uint64_t x[GF_MAX_DEGREE], root[GF_MAX_DEGREE], prod[GF_MAX_DEGREE];
	uint64_t q, k0, k;
	int d = 0, i;

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

	q = gf_field_order(p, n, GF_MAX_ORDER);

	if (q == 0) {
		return GF_ERR_NOT_FOUND;
	}

	/* alpha = x, a root of f. */
	memset(x, 0, n * sizeof(uint64_t));
	x[1] = 1;

	/* m(x) = product of (x - alpha^(e p^j)) over the distinct conjugates. */
	memset(m, 0, (n + 1) * n * sizeof(uint64_t));
	m[0] = 1;
	k0 = (uint64_t) e % (q - 1);
	k = k0;

	do {
		gf_elem_pow(p, n, f, x, k, root);

		/* m = m * (x - root) */
		memcpy(m + (d + 1) * n, m + d * n, n * sizeof(uint64_t));

		for (i = d; i > 0; i--) {
			gf_elem_mul(p, n, f, root, m + i * n, prod);
			gf_elem_sub(p, n, m + (i - 1) * n, prod, m + i * n);
		}

		gf_elem_mul(p, n, f, root, m, prod);
		memset(m, 0, n * sizeof(uint64_t));
		gf_elem_sub(p, n, m, prod, m);
		d++;
		k = (k * p) % (q - 1);
	} while (k != k0 && d < n);

	/* The coefficients lie in the prime subfield: only their constant term
	 * can be non-zero. */
	for (i = 0; i <= d; i++) {
		out[i] = m[i * n];
	}

	return GF_OK;
}

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
