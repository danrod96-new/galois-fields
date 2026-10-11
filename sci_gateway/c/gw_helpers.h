/*
 * gw_helpers.h - Argument checking and conversion shared by the
 * galois_fields gateways (Scilab 6+ API).
 *
 * Copyright (C) Daniel E Rodriguez - UIS.
 */
#ifndef GW_HELPERS_H
#define GW_HELPERS_H

#include <stdint.h>
#include "api_scilab.h"
#include "galois_fields.h"

/* A field as stored in Scilab: list(p, n, primitive_polynomial, elements). */
typedef struct
{
	uint64_t p;
	int n;
	uint64_t q;
	uint64_t f[GF_MAX_DEGREE + 1];
} gw_field;

/* Each helper reports its own error with Scierror and returns 1 on error. */
int gw_check_nin(const char *fname, int nin, int min, int max);
int gw_get_real(scilabEnv env, scilabVar var, const char* fname, int pos, double **data, int *rows, int *cols);
int gw_get_scalar_uint(scilabEnv env, scilabVar var, const char *fname, int pos, uint64_t *value);
int gw_get_scalar_int(scilabEnv env, scilabVar var, const char *fname, int pos, int64_t *value);
int gw_get_prime(scilabEnv env, scilabVar var, const char* fname, int pos, uint64_t* p);
int gw_get_degree(scilabEnv env, scilabVar var, const char* fname, int pos, int *n);
int gw_check_order(const char *fname, uint64_t p, int n, uint64_t *q);

/* Row vector of coefficients in 0 .. p - 1 (malloc'ed, caller frees). */
int gw_get_coeffs(scilabEnv env, scilabVar var, const char *fname, int pos, uint64_t p, uint64_t **coeffs, int* len);

/* Field list created by gf_create_field. */
int gw_get_field(scilabEnv env, scilabVar var, const char* fname, int pos, gw_field* field);

/* Elements of the field: for n = 1 any matrix of integers (count = its size),
 * for n > 1 a matrix with n columns, one element per row (count = rows).
 * data (malloc'ed) holds count * n coefficients, element after element. */
int gw_get_elements(scilabEnv env, scilabVar var, const char *fname, int pos, const gw_field *field, uint64_t **data, int *count, int *rows, int *cols);

/* Builds the Scilab value for count elements (n = 1: rows x cols matrix,
 * n > 1: count x n matrix). */
scilabVar gw_create_elements(scilabEnv env, const gw_field *field, const uint64_t *data, int count, int rows, int cols);

/* Row vector of len values. */
scilabVar gw_create_row(scilabEnv env, const uint64_t *values, int len);

/* Polynomial trimmed to its degree (the zero polynomial becomes 0). */
scilabVar gw_create_poly(scilabEnv env, const uint64_t *coeffs, int len);

#endif /* GW_HELPERS_H */
