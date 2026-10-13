/*
 * sci_gf_field.c - Gateways that create fields and work with primitive and
 * minimal polynomials.
 *
 * Copyright (C) Daniel E Rodriguez - UIS.
 */
#include <stdlib.h>
#include <string.h>
#include "api_scilab.h"
#include "Scierror.h"
#include "localization.h"
#include "gw_helpers.h"

/* field = gf_create_field(p, n [, primitive_polynomial]) */
int sci_gf_create_field(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_create_field";
	gw_field field;
	uint64_t* coeffs = NULL;
	int len = 0;
	scilabVar list, elements;
	double* table = NULL;

	if (gw_check_nin(fname, nin, 2, 3) || gw_get_prime(env, in[0], fname, 1, &field.p) || gw_get_degree(env, in[1], fname, 2, &field.n) || gw_check_order(fname, field.p, field.n, &field.q)) {
		return 1;
	}

	if (nin == 3) {
		if (gw_get_coeffs(env, in[2], fname, 3, field.p, &coeffs, &len)) {
			return 1;
		}

		if (len != field.n + 1 || coeffs[field.n] == 0) {
			free(coeffs);
			Scierror(999, _("%s: Wrong size for input argument #%d: A polynomial of degree %d expected (%d coefficients, constant term first).\n"), fname, 3, field.n, field.n + 1);
			return 1;
		}

		memcpy(field.f, coeffs, len * sizeof(uint64_t));
		free(coeffs);
		
		if (!gf_is_primitive_polynomial(field.f, field.n, field.p)) {
			Scierror(999, _("%s: Wrong value for input argument #%d: A primitive polynomial over GF(%llu) expected.\n"), fname, 3, (unsigned long long) field.p);
			return 1;
		}
	} else if (gf_default_primitive_polynomial(field.p, field.n, field.f) != GF_OK) {
		Scierror(999, _("%s: No primitive polynomial of degree %d found over GF(%llu).\n"), fname, field.n, (unsigned long long) field.p);
		return 1;
	}

	if (field.n == 1) {
		elements = scilab_createDoubleMatrix2d(env, 1, (int) field.q, 0);
	} else {
		elements = scilab_createDoubleMatrix2d(env, (int) field.q, field.n, 0);
	}

	scilab_getDoubleArray(env, elements, &table);
	gf_field_elements(field.p, field.n, field.f, table);

	list = scilab_createList(env);

	scilab_appendToList(env, list, scilab_createDouble(env, (double) field.p));
	scilab_appendToList(env, list, scilab_createDouble(env, (double) field.n));
	scilab_appendToList(env, list, gw_create_row(env, field.f, field.n + 1));
	scilab_appendToList(env, list, elements);

	out[0] = list;

	return 0;
}

/* q = gf_order(field) */
int sci_gf_order(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_order";
	gw_field field;

	if (gw_check_nin(fname, nin, 1, 1) || gw_get_field(env, in[0], fname, 1, &field)) {
		return 1;
	}

	out[0] = scilab_createDouble(env, (double)field.q);
	
	return 0;
}

/* f = gf_default_primitive_polynomial(p, n) */
int sci_gf_default_primitive_polynomial(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_default_primitive_polynomial";
	uint64_t p = 0, q = 0;
	uint64_t f[GF_MAX_DEGREE + 1];
	int n = 0;

	if (gw_check_nin(fname, nin, 2, 2) || gw_get_prime(env, in[0], fname, 1, &p) || gw_get_degree(env, in[1], fname, 2, &n) || gw_check_order(fname, p, n, &q)) {
		return 1;
	}

	if (gf_default_primitive_polynomial(p, n, f) != GF_OK) {
		Scierror(999, _("%s: No primitive polynomial of degree %d found over GF(%llu).\n"), fname, n, (unsigned long long) p);
		return 1;
	}

	out[0] = gw_create_row(env, f, n + 1);
	
	return 0;
}

/* tf = gf_is_primitive_polynomial(f, p, n) */
int sci_gf_is_primitive_polynomial(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
	static const char fname[] = "gf_is_primitive_polynomial";
	uint64_t p = 0, q = 0;
	uint64_t* f = NULL;
	int n = 0, len = 0, result;

	if (gw_check_nin(fname, nin, 3, 3) || gw_get_prime(env, in[1], fname, 2, &p) || gw_get_degree(env, in[2], fname, 3, &n) || gw_check_order(fname, p, n, &q) || gw_get_coeffs(env, in[0], fname, 1, p, &f, &len)) {
		return 1;
	}

	if (len != n + 1 || f[n] == 0) {
		free(f);
		Scierror(999, _("%s: Wrong size for input argument #%d: A polynomial of degree %d expected (%d coefficients, constant term first).\n"), fname, 1, n, n + 1);
		return 1;
	}
	
	result = gf_is_primitive_polynomial(f, n, p);
	free(f);
	out[0] = scilab_createBoolean(env, result);
	
	return 0;
}

/* m = gf_minimal_polynomial(p, n, e) or gf_minimal_polynomial(field, e)
 * e is an element of GF(p) for n = 1, or the exponent of alpha^e for n > 1
 * (-1 for the zero element). */
int sci_gf_minimal_polynomial(scilabEnv env, int nin, scilabVar *in, int nopt, scilabOpt *opt, int nout, scilabVar *out)
{
	static const char fname[] = "gf_minimal_polynomial";
	gw_field field;
	int64_t e = 0;
	int epos;
	uint64_t result[GF_MAX_DEGREE + 1];
	int status;

	if (gw_check_nin(fname, nin, 2, 3)) {
		return 1;
	}

	if (scilab_isList(env, in[0])) {
		if (nin != 2) {
			Scierror(77, _("%s: Wrong number of input argument(s): %d expected.\n"), fname, 2);
			return 1;
		}

		if (gw_get_field(env, in[0], fname, 1, &field)) {
			return 1;
		}

		epos = 2;
	} else {
		if (nin != 3 || gw_get_prime(env, in[0], fname, 1, &field.p) || gw_get_degree(env, in[1], fname, 2, &field.n) || gw_check_order(fname, field.p, field.n, &field.q)) {
			if (nin != 3) {
				Scierror(77, _("%s: Wrong number of input argument(s): %d expected.\n"), fname, 3);
			}
			
			return 1;
		}

		gf_default_primitive_polynomial(field.p, field.n, field.f);
		epos = 3;
	}

	if (gw_get_scalar_int(env, in[epos - 1], fname, epos, &e)) {
		return 1;
	}

	if (field.n == 1 && (e < 0 || (uint64_t) e >= field.p)) {
		Scierror(999, _("%s: Wrong value for input argument #%d: An element in [0, %llu] expected.\n"), fname, epos, (unsigned long long) (field.p - 1));
		return 1;
	}

	if (field.n > 1 && (e < -1 || (uint64_t) (e + 1) >= field.q)) {
		Scierror(999, _("%s: Wrong value for input argument #%d: An exponent in [-1, %llu] expected (-1 is the zero element).\n"), fname, epos, (unsigned long long) (field.q - 2));
		return 1;
	}

	status = gf_minimal_polynomial(field.p, field.n, field.f, e, result);

	if (status != GF_OK) {
		Scierror(999, _("%s: Memory allocation error.\n"), fname);
		return 1;
	}

	out[0] = gw_create_row(env, result, field.n + 1);
	
	return 0;
}
