// Unit tests: polynomials over GF(p) and integer helpers.
// <-- NO CHECK REF -->
// <-- CLI SHELL MODE -->

assert_checkequal(gf_poly_add([1 2], [2 2 1], 3), [0 1 1]);
assert_checkequal(gf_poly_add([1 1], [2 2], 3), 0);                  // zero polynomial
assert_checkequal(gf_poly_subtract([1 2], [2 2 1], 3), [2 0 2]);
assert_checkequal(gf_poly_multiply([1 1], [1 1], 2), [1 0 1]);       // (x + 1)^2 = x^2 + 1
assert_checkequal(gf_poly_multiply([1 2], [1 1], 3), [1 0 2]);       // (2x + 1)(x + 1) = 2x^2 + 1
assert_checkequal(gf_poly_multiply([0], [1 1], 3), 0);
assert_checkequal(gf_poly_power([1 1], 2, 2), [1 0 1]);
assert_checkequal(gf_poly_power([1 1], 8, 2), [1 zeros(1, 7) 1]);    // Frobenius
assert_checkequal(gf_poly_power([1 1], 0, 2), 1);
assert_checkequal(gf_poly_power([0], 3, 2), 0);

// Division: a = q * b + r with deg r < deg b
a = [1 0 0 0 0 0 1];   // x^6 + 1
b = [1 1 0 1];         // x^3 + x + 1
[r, q] = gf_poly_divide(a, b, 2);
assert_checkequal(q, [1 1 0 1]);
assert_checkequal(r, [0 0 1]);
assert_checkequal(gf_poly_add(gf_poly_multiply(q, b, 2), r, 2), a);

a = [4 0 3 1 2];       // over GF(5), non-monic divisor
b = [1 3 2];
[r, q] = gf_poly_divide(a, b, 5);
assert_checktrue(size(r, "*") < 3);
assert_checkequal(gf_poly_add(gf_poly_multiply(q, b, 5), r, 5), a);

[r, q] = gf_poly_divide([1 1], [1 0 1], 2);   // degree of a below degree of b
assert_checkequal(q, 0);
assert_checkequal(r, [1 1]);

// Integer helpers
assert_checkequal(gf_is_prime([0 1 2 3 4 97 91 1048573]), [%f %f %t %t %f %t %f %t]);
assert_checkequal(gf_mod([-3 10 0; 7 -14 13], 7), [4 3 0; 0 0 6]);

// Errors
assert_checkerror("gf_poly_divide([1 1], 0, 2)", [], 999);
assert_checkerror("gf_poly_add([1 2], [1 1], 2)", [], 999);   // coefficient >= p
assert_checkerror("gf_poly_add([1 1], [1 1], 4)", [], 999);   // 4 is not prime
assert_checkerror("gf_poly_add([1; 1], [1 1], 2)", [], 999);  // column vector
assert_checkerror("gf_poly_power([0], 0, 2)", [], 999);
assert_checkerror("gf_is_prime(-1)", [], 999);
assert_checkerror("gf_mod(1.5, 7)", [], 999);
assert_checkerror("gf_mod(3, 0)", [], 999);
