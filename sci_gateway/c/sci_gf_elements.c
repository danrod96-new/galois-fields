/*
 * sci_gf_elements.c - Gateways for arithmetic between field elements.
 *
 * Elements of GF(p) are integers; elements of GF(p^n), n > 1, are rows of n
 * coefficients (constant term first). Several elements can be given at once
 * (a matrix for n = 1, one element per row for n > 1); a single element is
 * combined with each of the others.
 *
 * Copyright (C) Daniel E Rodriguez - UIS.
 */
#include <stdlib.h>
#include "api_scilab.h"
#include "Scierror.h"
#include "localization.h"
#include "gw_helpers.h"

enum gw_op { GW_ADD, GW_SUB, GW_MUL, GW_DIV };

static int binary_op(scilabEnv env, int nin, scilabVar* in, scilabVar* out, const char* fname, enum gw_op op)
{
    gw_field field;
    uint64_t *a = NULL, *b = NULL, *r = NULL;
    int ca = 0, cb = 0, ra = 0, rb = 0, cola = 0, colb = 0;
    int count, rows, cols, i, n, status = 0;

    if (gw_check_nin(fname, nin, 3, 3)
            || gw_get_field(env, in[0], fname, 1, &field)
            || gw_get_elements(env, in[1], fname, 2, &field, &a, &ca, &ra, &cola))
    {
        return 1;
    }
    if (gw_get_elements(env, in[2], fname, 3, &field, &b, &cb, &rb, &colb))
    {
        free(a);
        return 1;
    }
    if (ca != cb && ca != 1 && cb != 1)
    {
        free(a);
        free(b);
        Scierror(999, _("%s: Wrong size for input arguments #%d and #%d: Same number of elements expected, or a single element.\n"),
                 fname, 2, 3);
        return 1;
    }

    n = field.n;
    count = ca >= cb ? ca : cb;
    rows = ca >= cb ? ra : rb;
    cols = ca >= cb ? cola : colb;
    r = (uint64_t*)malloc((size_t)count * n * sizeof(uint64_t));
    if (r == NULL)
    {
        free(a);
        free(b);
        Scierror(999, _("%s: Memory allocation error.\n"), fname);
        return 1;
    }

    for (i = 0; i < count && status == 0; i++)
    {
        const uint64_t* x = a + (ca == 1 ? 0 : i) * n;
        const uint64_t* y = b + (cb == 1 ? 0 : i) * n;
        uint64_t* z = r + i * n;
        switch (op)
        {
            case GW_ADD:
                gf_elem_add(field.p, n, x, y, z);
                break;
            case GW_SUB:
                gf_elem_sub(field.p, n, x, y, z);
                break;
            case GW_MUL:
                gf_elem_mul(field.p, n, field.f, x, y, z);
                break;
            case GW_DIV:
                status = gf_elem_div(field.p, n, field.f, x, y, z);
                break;
        }
    }

    if (status != 0)
    {
        Scierror(999, _("%s: Division by zero.\n"), fname);
    }
    else
    {
        out[0] = gw_create_elements(env, &field, r, count, rows, cols);
    }
    free(a);
    free(b);
    free(r);
    return status != 0;
}

/* c = gf_add(field, a, b) */
int sci_gf_add(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
    return binary_op(env, nin, in, out, "gf_add", GW_ADD);
}

/* c = gf_subtract(field, a, b) */
int sci_gf_subtract(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
    return binary_op(env, nin, in, out, "gf_subtract", GW_SUB);
}

/* c = gf_multiply(field, a, b) */
int sci_gf_multiply(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
    return binary_op(env, nin, in, out, "gf_multiply", GW_MUL);
}

/* c = gf_divide(field, a, b) */
int sci_gf_divide(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
    return binary_op(env, nin, in, out, "gf_divide", GW_DIV);
}

/* c = gf_power(field, a, e), e a non-negative integer */
int sci_gf_power(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
    static const char fname[] = "gf_power";
    gw_field field;
    uint64_t *a = NULL, *r = NULL, e = 0;
    int count = 0, rows = 0, cols = 0, i, status = 0;

    if (gw_check_nin(fname, nin, 3, 3)
            || gw_get_field(env, in[0], fname, 1, &field)
            || gw_get_scalar_uint(env, in[2], fname, 3, &e)
            || gw_get_elements(env, in[1], fname, 2, &field, &a, &count, &rows, &cols))
    {
        return 1;
    }
    r = (uint64_t*)malloc((size_t)count * field.n * sizeof(uint64_t));
    if (r == NULL)
    {
        free(a);
        Scierror(999, _("%s: Memory allocation error.\n"), fname);
        return 1;
    }
    for (i = 0; i < count && status == 0; i++)
    {
        status = gf_elem_pow(field.p, field.n, field.f, a + i * field.n, e, r + i * field.n);
    }
    if (status != 0)
    {
        Scierror(999, _("%s: 0^0 is undefined.\n"), fname);
    }
    else
    {
        out[0] = gw_create_elements(env, &field, r, count, rows, cols);
    }
    free(a);
    free(r);
    return status != 0;
}
