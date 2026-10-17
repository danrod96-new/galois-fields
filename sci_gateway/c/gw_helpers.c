/*
 * gw_helpers.c - Argument checking and conversion shared by the
 * galois_fields gateways (Scilab 6+ API).
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

#define GW_MAX_EXACT 9007199254740992.0 /* 2^53 */

int gw_check_nin(const char *fname, int nin, int min, int max)
{
	if (nin < min || nin > max)
	{
		if (min == max) {
			Scierror(77, _("%s: Wrong number of input argument(s): %d expected.\n"), fname, min);
		} else {
			Scierror(77, _("%s: Wrong number of input argument(s): %d to %d expected.\n"), fname, min, max);
		}
		return 1;
	}

	return 0;
}

int gw_get_real(scilabEnv env, scilabVar var, const char* fname, int pos, double **data, int *rows, int *cols)
{
	if (scilab_isDouble(env, var) == 0 || scilab_isMatrix2d(env, var) == 0 || scilab_isComplex(env, var) == 1) {
		Scierror(999, _("%s: Wrong type for input argument #%d: A real matrix expected.\n"), fname, pos);
		return 1;
	}

	scilab_getDim2d(env, var, rows, cols);
	scilab_getDoubleArray(env, var, data);

	return 0;
}

static int is_integer(double v)
{
  return isfinite(v) && floor(v) == v && fabs(v) < GW_MAX_EXACT;
}

int gw_get_scalar_int(scilabEnv env, scilabVar var, const char *fname, int pos, int64_t *value)
{
	double *data = NULL;
	int rows = 0, cols = 0;
	
	if (gw_get_real(env, var, fname, pos, &data, &rows, &cols)) {
		return 1;
	}

	if (rows * cols != 1 || !is_integer(data[0])) {
		Scierror(999, _("%s: Wrong value for input argument #%d: An integer scalar expected.\n"), fname, pos);
		return 1;
	}

	*value = (int64_t) data[0];
	
	return 0;
}

int gw_get_scalar_uint(scilabEnv env, scilabVar var, const char *fname, int pos, uint64_t *value)
{
	int64_t v = 0;

	if (gw_get_scalar_int(env, var, fname, pos, &v)) {
		return 1;
	}

	if (v < 0) {
		Scierror(999, _("%s: Wrong value for input argument #%d: A non-negative integer expected.\n"), fname, pos);
		return 1;
	}

	*value = (uint64_t) v;
	
	return 0;
}

int gw_get_prime(scilabEnv env, scilabVar var, const char* fname, int pos, uint64_t *p)
{
	if (gw_get_scalar_uint(env, var, fname, pos, p)) {
		return 1;
	}

	if (!gf_is_prime(*p)) {
		Scierror(999, _("%s: Wrong value for input argument #%d: A prime number expected.\n"), fname, pos);
		return 1;
	}

	if (*p > GF_MAX_ORDER) {
		Scierror(999, _("%s: Wrong value for input argument #%d: The characteristic must not exceed %llu.\n"), fname, pos, (unsigned long long) GF_MAX_ORDER);
		return 1;
	}

	return 0;
}

int gw_get_degree(scilabEnv env, scilabVar var, const char *fname, int pos, int *n)
{
	uint64_t v = 0;

	if (gw_get_scalar_uint(env, var, fname, pos, &v)) {
		return 1;
	}

	if (v < 1 || v > GF_MAX_DEGREE) {
		Scierror(999, _("%s: Wrong value for input argument #%d: An integer in [%d, %d] expected.\n"), fname, pos, 1, GF_MAX_DEGREE);
		return 1;
	}

	*n = (int) v;
	
	return 0;
}

int gw_check_order(const char *fname, uint64_t p, int n, uint64_t *q)
{
	*q = gf_field_order(p, n, GF_MAX_ORDER);

	if (*q == 0) {
		Scierror(999, _("%s: The field order p^n must not exceed %llu.\n"), fname, (unsigned long long) GF_MAX_ORDER);
		return 1;
	}

	return 0;
}

static int read_coeff(double v, uint64_t p, uint64_t *out)
{
	if (!is_integer(v) || v < 0 || (uint64_t)v >= p) {
		return 1;
	}

	*out = (uint64_t) v;
	
	return 0;
}

int gw_get_coeffs(scilabEnv env, scilabVar var, const char *fname, int pos, uint64_t p, uint64_t **coeffs, int *len)
{
	double *data = NULL;
	int rows = 0, cols = 0, i;

	if (gw_get_real(env, var, fname, pos, &data, &rows, &cols)) {
		return 1;
	}

	if (rows != 1 || cols < 1) {
		Scierror(999, _("%s: Wrong size for input argument #%d: A row vector of coefficients expected.\n"), fname, pos);
		return 1;
	}

	*coeffs = (uint64_t *) malloc(cols * sizeof(uint64_t));

	if (*coeffs == NULL) {
		Scierror(999, _("%s: Memory allocation error.\n"), fname);
		return 1;
	}

	for (i = 0; i < cols; i++) {
		if (read_coeff(data[i], p, &(*coeffs)[i])) {
			free(*coeffs);
			*coeffs = NULL;
			Scierror(999, _("%s: Wrong value for input argument #%d: Integer coefficients in [0, %llu] expected.\n"), fname, pos, (unsigned long long) (p - 1));
			return 1;
		}
	}

	*len = cols;

	return 0;
}

static int field_error(const char *fname, int pos)
{
  Scierror(999, _("%s: Wrong type for input argument #%d: A field created by gf_create_field expected.\n"), fname, pos);
  
	return 1;
}

int gw_get_field(scilabEnv env, scilabVar var, const char *fname, int pos, gw_field *field)
{
	scilabVar item;
	double *data = NULL;
	int rows = 0, cols = 0, i;

	if (scilab_isList(env, var) == 0 || scilab_getSize(env, var) != 4) {
		return field_error(fname, pos);
	}

	item = scilab_getListItem(env, var, 0);

	if (scilab_isDouble(env, item) == 0 || scilab_getSize(env, item) != 1) {
		return field_error(fname, pos);
	}

	scilab_getDoubleArray(env, item, &data);

	if (!is_integer(data[0]) || !gf_is_prime((uint64_t) data[0]) || data[0] > (double) GF_MAX_ORDER) {
		return field_error(fname, pos);
	}

	field->p = (uint64_t) data[0];

	item = scilab_getListItem(env, var, 1);
	
	if (scilab_isDouble(env, item) == 0 || scilab_getSize(env, item) != 1) {
		return field_error(fname, pos);
	}

	scilab_getDoubleArray(env, item, &data);
	
	if (!is_integer(data[0]) || data[0] < 1 || data[0] > GF_MAX_DEGREE) {
		return field_error(fname, pos);
	}

	field->n = (int) data[0];
	field->q = gf_field_order(field->p, field->n, GF_MAX_ORDER);
	
	if (field->q == 0) {
		return field_error(fname, pos);
	}

	item = scilab_getListItem(env, var, 2);
	
	if (scilab_isDouble(env, item) == 0 || scilab_isComplex(env, item) == 1) {
		return field_error(fname, pos);
	}

	scilab_getDim2d(env, item, &rows, &cols);
	
	if (rows != 1 || cols != field->n + 1) {
		return field_error(fname, pos);
	}

	scilab_getDoubleArray(env, item, &data);
	
	for (i = 0; i <= field->n; i++) {
		if (read_coeff(data[i], field->p, &field->f[i])) {
			return field_error(fname, pos);
		}
	}
	
	if (field->f[field->n] == 0) {
		return field_error(fname, pos);
	}

	item = scilab_getListItem(env, var, 3);

	if (scilab_isDouble(env, item) == 0) {
		return field_error(fname, pos);
	}

	scilab_getDim2d(env, item, &rows, &cols);
	
	if ((uint64_t) rows * (uint64_t) cols != field->q * (uint64_t) field->n) {
		return field_error(fname, pos);
	}

	return 0;
}

int gw_get_elements(scilabEnv env, scilabVar var, const char *fname, int pos, const gw_field *field, uint64_t **data, int *count, int *rows, int *cols)
{
  double *values = NULL;
  int r, c;
  int n = field->n;

  if (gw_get_real(env, var, fname, pos, &values, rows, cols)) {
    return 1;
  }

  if (*rows * *cols == 0) {
    Scierror(999, _("%s: Wrong size for input argument #%d: A non-empty matrix expected.\n"), fname, pos);
    return 1;
  }

  if (n == 1) {
    *count = *rows * *cols;
  } else {
    if (*cols != n) {
      Scierror(999, _("%s: Wrong size for input argument #%d: %d columns expected (one element per row).\n"), fname, pos, n);
      return 1;
    }
    
    *count = *rows;
  }

  *data = (uint64_t *) malloc((size_t) *count * n * sizeof(uint64_t));
  
  if (*data == NULL) {
    Scierror(999, _("%s: Memory allocation error.\n"), fname);
    return 1;
  }

  for (r = 0; r < *count; r++) {
    for (c = 0; c < n; c++) {
      /* Scilab stores matrices column by column. */
      double v = n == 1 ? values[r] : values[c * *rows + r];

      if (read_coeff(v, field->p, &(*data)[r * n + c])) {
        free(*data);
        *data = NULL;
        Scierror(999, _("%s: Wrong value for input argument #%d: Integers in [0, %llu] expected.\n"), fname, pos, (unsigned long long) (field->p - 1));
        return 1;
      }
    }
  }
  
  return 0;
}

scilabVar gw_create_elements(scilabEnv env, const gw_field *field, const uint64_t *data, int count, int rows, int cols)
{
	scilabVar out;
	double *values = NULL;
	int r, c, n = field->n;

	if (n == 1) {
		out = scilab_createDoubleMatrix2d(env, rows, cols, 0);
		scilab_getDoubleArray(env, out, &values);

		for (r = 0; r < count; r++) {
			values[r] = (double)data[r];
		}

		return out;
	}

	out = scilab_createDoubleMatrix2d(env, count, n, 0);
	scilab_getDoubleArray(env, out, &values);
	
	for (r = 0; r < count; r++) {
		for (c = 0; c < n; c++) {
			values[c * count + r] = (double)data[r * n + c];
		}
	}

	return out;
}

scilabVar gw_create_row(scilabEnv env, const uint64_t *values, int len)
{
	scilabVar out = scilab_createDoubleMatrix2d(env, 1, len, 0);
	double *data = NULL;
	int i;

	scilab_getDoubleArray(env, out, &data);
	
	for (i = 0; i < len; i++) {
		data[i] = (double) values[i];
	}

	return out;
}

scilabVar gw_create_poly(scilabEnv env, const uint64_t *coeffs, int len)
{
	int degree = gf_poly_degree(coeffs, len);

	if (degree < 0) {
		return scilab_createDouble(env, 0);
	}

	return gw_create_row(env, coeffs, degree + 1);
}
