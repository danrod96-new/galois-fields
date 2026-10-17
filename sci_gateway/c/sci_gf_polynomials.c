/*
 * sci_gf_polynomials.c - Gateways for polynomials over GF(p) and for the
 * integer helpers.
 *
 * Polynomials are row vectors of coefficients in 0 .. p - 1, constant term
 * first. Results are trimmed to their degree; the zero polynomial is 0.
 *
 * Copyright (C) Daniel E Rodriguez - UIS.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "api_scilab.h"
#include "Scierror.h"
#include "localization.h"
#include "gw_helpers.h"

enum gw_poly_op { GW_PADD, GW_PSUB, GW_PMUL };

static int read_two_polys(scilabEnv env, int nin, scilabVar *in, const char *fname, uint64_t *p, uint64_t **a, int *la, uint64_t **b, int *lb)
{
	if (gw_check_nin(fname, nin, 3, 3) || gw_get_prime(env, in[2], fname, 3, p) || gw_get_coeffs(env, in[0], fname, 1, *p, a, la)) {
		return 1;
	}

	if (gw_get_coeffs(env, in[1], fname, 2, *p, b, lb)) {
		free(*a);
		return 1;
	}

	return 0;
}

static int poly_op(scilabEnv env, int nin, scilabVar *in, scilabVar *out, const char *fname, enum gw_poly_op op)
{
	uint64_t p = 0, *a = NULL, *b = NULL, *r = NULL;
	int la = 0, lb = 0, lr;

	if (read_two_polys(env, nin, in, fname, &p, &a, &la, &b, &lb)) {
		return 1;
	}

	lr = op == GW_PMUL ? la + lb - 1 : (la > lb ? la : lb);
	r = (uint64_t*) malloc(lr * sizeof(uint64_t));
	
	if (r == NULL) {
		free(a);
		free(b);
		Scierror(999, _("%s: Memory allocation error.\n"), fname);
		return 1;
	}

	switch (op) {
		case GW_PADD:
			gf_poly_add(p, a, la, b, lb, r);
			break;
		case GW_PSUB:
			gf_poly_sub(p, a, la, b, lb, r);
			break;
		case GW_PMUL:
			gf_poly_mul(p, a, la, b, lb, r);
			break;
	}

	out[0] = gw_create_poly(env, r, lr);

	free(a);
	free(b);
	free(r);

	return 0;
}

/* c = gf_poly_add(a, b, p) */
int sci_gf_poly_add(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
  return poly_op(env, nin, in, out, "gf_poly_add", GW_PADD);
}

/* c = gf_poly_subtract(a, b, p) */
int sci_gf_poly_subtract(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
  return poly_op(env, nin, in, out, "gf_poly_subtract", GW_PSUB);
}

/* c = gf_poly_multiply(a, b, p) */
int sci_gf_poly_multiply(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
  return poly_op(env, nin, in, out, "gf_poly_multiply", GW_PMUL);
}

/* [r, q] = gf_poly_divide(a, b, p): a = q * b + r, deg r < deg b */
int sci_gf_poly_divide(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_poly_divide";
	uint64_t p = 0, *a = NULL, *b = NULL, *q = NULL, *r = NULL;
	int la = 0, lb = 0;

	if (nout > 2) {
		Scierror(78, _("%s: Wrong number of output argument(s): %d to %d expected.\n"), fname, 1, 2);
		return 1;
	}

	if (read_two_polys(env, nin, in, fname, &p, &a, &la, &b, &lb)) {
		return 1;
	}

	q = (uint64_t *) malloc(la * sizeof(uint64_t));
	r = (uint64_t *) malloc(la * sizeof(uint64_t));
	
	if (q == NULL || r == NULL) {
		free(a);
		free(b);
		free(q);
		free(r);
		Scierror(999, _("%s: Memory allocation error.\n"), fname);
		return 1;
	}

	if (gf_poly_divmod(p, a, la, b, lb, q, r) != GF_OK) {
		free(a);
		free(b);
		free(q);
		free(r);
		Scierror(999, _("%s: Division by the zero polynomial.\n"), fname);
		return 1;
	}

	out[0] = gw_create_poly(env, r, la);

	if (nout == 2) {
		out[1] = gw_create_poly(env, q, la);
	}

	free(a);
	free(b);
	free(q);
	free(r);

	return 0;
}

/* c = gf_poly_power(a, e, p) */
int sci_gf_poly_power(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_poly_power";
	uint64_t p = 0, e = 0, *a = NULL, *result = NULL, *base = NULL, *tmp = NULL;
	int la = 0, lr = 1, lbase, deg;

	if (gw_check_nin(fname, nin, 3, 3) || gw_get_prime(env, in[2], fname, 3, &p) || gw_get_scalar_uint(env, in[1], fname, 2, &e) || gw_get_coeffs(env, in[0], fname, 1, p, &a, &la)) {
		return 1;
	}

	deg = gf_poly_degree(a, la);

	if (deg < 0) {
		free(a);

		if (e == 0) {
			Scierror(999, _("%s: 0^0 is undefined.\n"), fname);
			return 1;
		}

		out[0] = scilab_createDouble(env, 0);
		return 0;
	}

	if ((uint64_t) deg * e > 100000) {
		free(a);
		Scierror(999, _("%s: The result would have degree %llu; at most %d is supported.\n"), fname, (unsigned long long)((uint64_t)deg * e), 100000);
		return 1;
	}

	lbase = deg + 1;
	result = (uint64_t *) calloc((size_t) deg * e + 1, sizeof(uint64_t));
	base = (uint64_t *) calloc((size_t) deg * e + 1, sizeof(uint64_t));
	tmp = (uint64_t *) calloc((size_t) deg * e + 1, sizeof(uint64_t));

	if (result == NULL || base == NULL || tmp == NULL) {
		free(a);
		free(result);
		free(base);
		free(tmp);
		Scierror(999, _("%s: Memory allocation error.\n"), fname);
		return 1;
	}

	result[0] = 1;
	memcpy(base, a, lbase * sizeof(uint64_t));

	while (e > 0) {
		if (e & 1) {
			gf_poly_mul(p, result, lr, base, lbase, tmp);
			lr = lr + lbase - 1;
			memcpy(result, tmp, lr * sizeof(uint64_t));
		}

		e >>= 1;

		if (e > 0) {
			gf_poly_mul(p, base, lbase, base, lbase, tmp);
			lbase = 2 * lbase - 1;
			memcpy(base, tmp, lbase * sizeof(uint64_t));
		}
	}

	out[0] = gw_create_poly(env, result, lr);
	
	free(a);
	free(result);
	free(base);
	free(tmp);
	
	return 0;
}

/* tf = gf_is_prime(x), element-wise on a matrix of non-negative integers */
int sci_gf_is_prime(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_is_prime";
	double *data = NULL;
	int *res = NULL;
	int rows = 0, cols = 0, i;

	if (gw_check_nin(fname, nin, 1, 1) || gw_get_real(env, in[0], fname, 1, &data, &rows, &cols)) {
		return 1;
	}

	for (i = 0; i < rows * cols; i++) {
		if (!(data[i] >= 0 && floor(data[i]) == data[i] && data[i] < 9007199254740992.0)) {
			Scierror(999, _("%s: Wrong value for input argument #%d: Non-negative integers expected.\n"), fname, 1);
      return 1;
		}
	}

	out[0] = scilab_createBooleanMatrix2d(env, rows, cols);
	scilab_getBooleanArray(env, out[0], &res);

	for (i = 0; i < rows * cols; i++) {
		res[i] = gf_is_prime((uint64_t)data[i]);
	}

	return 0;
}

/* r = gf_mod(x, m): non-negative remainder of the integers x modulo m > 0 */
int sci_gf_mod(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_mod";
	double *data = NULL, *res = NULL;
	uint64_t m = 0;
	int rows = 0, cols = 0, i;

	if (gw_check_nin(fname, nin, 2, 2) || gw_get_real(env, in[0], fname, 1, &data, &rows, &cols) || gw_get_scalar_uint(env, in[1], fname, 2, &m)) {
		return 1;
	}

	if (m == 0){
		Scierror(999, _("%s: Wrong value for input argument #%d: A positive integer expected.\n"), fname, 2);
    return 1;
	}

	for (i = 0; i < rows * cols; i++) {
		if (!(floor(data[i]) == data[i] && fabs(data[i]) < 9007199254740992.0)){
			Scierror(999, _("%s: Wrong value for input argument #%d: Integers expected.\n"), fname, 1);
      return 1;
		}
	}
	out[0] = scilab_createDoubleMatrix2d(env, rows, cols, 0);
	scilab_getDoubleArray(env, out[0], &res);

	for (i = 0; i < rows * cols; i++) {
		res[i] = (double) gf_mod((int64_t) data[i], m);
	}
	
	return 0;
}
