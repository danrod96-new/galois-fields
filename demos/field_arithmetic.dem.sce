// A short tour of GF(2^4) and GF(7).
// Copyright (C) Daniel E Rodriguez - UIS.
function demo_field_arithmetic()
    mode(-1);
    lines(0);

    mprintf("GF(2^4) with the primitive polynomial 1 + x + x^4\n");
    f = gf_create_field(2, 4);
    disp(gf_primitive_polynomial(f));

    E = f(4);
    mprintf("\nElements (row r is alpha^(r - 2); the first row is 0):\n");
    disp(E);

    a = E(6, :);   // alpha^4
    b = E(9, :);   // alpha^7
    mprintf("\na = alpha^4 = %s, b = alpha^7 = %s\n", sci2exp(a), sci2exp(b));
    mprintf("a + b = %s\n", sci2exp(gf_add(f, a, b)));
    mprintf("a * b = %s = alpha^%d\n", sci2exp(gf_multiply(f, a, b)), gf_to_exponential(f, gf_multiply(f, a, b)));
    mprintf("a / b = %s = alpha^%d\n", sci2exp(gf_divide(f, a, b)), gf_to_exponential(f, gf_divide(f, a, b)));
    mprintf("a^5   = %s\n", sci2exp(gf_power(f, a, 5)));

    mprintf("\nMinimal polynomial of alpha^5: %s\n", sci2exp(gf_minimal_polynomial(f, 5)));
    mprintf("Primitive polynomials of degree 4 over GF(2):\n");
    disp(gf_primitive_polynomials(2, 4));

    mprintf("\nMultiplication table of GF(7):\n");
    disp(gf_multiplication_table(gf_create_field(7, 1)));
endfunction

demo_field_arithmetic();
clear demo_field_arithmetic;
