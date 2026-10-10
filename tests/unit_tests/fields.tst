// Unit tests: creating fields, primitive and minimal polynomials.
// <-- NO CHECK REF -->
// <-- CLI SHELL MODE -->

// GF(2^4) with the default primitive polynomial x^4 + x + 1
f = gf_create_field(2, 4);
assert_checkequal(typeof(f), "list");
assert_checkequal(size(f), 4);
assert_checkequal(f(1), 2);
assert_checkequal(f(2), 4);
assert_checkequal(f(3), [1 1 0 0 1]);
E = f(4);
assert_checkequal(size(E), [16 4]);
assert_checkequal(E(1, :), [0 0 0 0]);      // zero
assert_checkequal(E(2, :), [1 0 0 0]);      // alpha^0
assert_checkequal(E(3, :), [0 1 0 0]);      // alpha^1
assert_checkequal(E(6, :), [1 1 0 0]);      // alpha^4 = 1 + alpha
assert_checkequal(size(unique(E, "r"), "r"), 16);
assert_checkequal(gf_order(f), 16);
assert_checkequal(gf_characteristic(f), 2);
assert_checkequal(gf_degree(f), 4);

// Explicit primitive polynomial (the one used by the Hamming demo)
g = gf_create_field(2, 4, [1 0 0 1 1]);
assert_checkequal(g(3), [1 0 0 1 1]);
assert_checkequal(size(unique(g(4), "r"), "r"), 16);

// Prime field GF(7): x - 3, 3 being the smallest primitive root
h = gf_create_field(7, 1);
assert_checkequal(h(3), [4 1]);
assert_checkequal(h(4), 0:6);
assert_checkequal(gf_order(h), 7);

// GF(3^2)
k = gf_create_field(3, 2);
assert_checkequal(k(3), [2 1 1]);
assert_checkequal(size(unique(k(4), "r"), "r"), 9);

// Default and primitive polynomials
assert_checkequal(gf_default_primitive_polynomial(2, 3), [1 1 0 1]);
assert_checkequal(gf_default_primitive_polynomial(2, 8), [1 0 1 1 1 0 0 0 1]);
assert_checkequal(gf_default_primitive_polynomial(5, 1), [3 1]);
assert_checktrue(gf_is_primitive_polynomial([1 0 0 1 1], 2, 4));
assert_checkfalse(gf_is_primitive_polynomial([1 1 1 1 1], 2, 4)); // irreducible, order 5
assert_checkfalse(gf_is_primitive_polynomial([1 0 1 0 1], 2, 4)); // reducible
assert_checkequal(gf_primitive_polynomials(2, 4), [1 1 0 0 1; 1 0 0 1 1]);
assert_checkequal(gf_primitive_polynomials(7, 1), [4 1; 2 1]);
assert_checkequal(size(gf_primitive_polynomials(5, 3), "r"), 20);  // phi(124) / 3
assert_checkequal(size(gf_primitive_polynomials(3, 4), "r"), 8);   // phi(80) / 4

// Minimal polynomials in GF(2^4) (x^4 + x + 1)
assert_checkequal(gf_minimal_polynomial(2, 4, 1), [1 1 0 0 1]);
assert_checkequal(gf_minimal_polynomial(2, 4, 0), [1 1 0 0 0]);   // of 1: x + 1
assert_checkequal(gf_minimal_polynomial(2, 4, 3), [1 1 1 1 1]);
assert_checkequal(gf_minimal_polynomial(2, 4, 5), [1 1 1 0 0]);   // degree 2
assert_checkequal(gf_minimal_polynomial(2, 4, -1), [0 1 0 0 0]);  // of 0: x
assert_checkequal(gf_minimal_polynomial(f, 5), [1 1 1 0 0]);
assert_checkequal(gf_minimal_polynomial(7, 1, 3), [4 1]);
// The minimal polynomial of alpha^e uses the field's own primitive polynomial
assert_checkequal(gf_minimal_polynomial(g, 1), [1 0 0 1 1]);

// Errors
assert_checkerror("gf_create_field(4, 2)", [], 999);         // 4 is not prime
assert_checkerror("gf_create_field(2, 0)", [], 999);
assert_checkerror("gf_create_field(2, 21)", [], 999);
assert_checkerror("gf_create_field(1048573, 2)", [], 999);   // p^n too large
assert_checkerror("gf_create_field(2, 4, [1 1 1 1 1])", [], 999);
assert_checkerror("gf_create_field(2, 4, [1 1 0 1])", [], 999);
assert_checkerror("gf_create_field(2, 4, [1 2 0 0 1])", [], 999);
assert_checkerror("gf_create_field(2.5, 4)", [], 999);
assert_checkerror("gf_create_field(2)", [], 77);
assert_checkerror("gf_order(list(1, 2))", [], 999);
assert_checkerror("gf_degree(42)", [], 10000);
assert_checkerror("gf_minimal_polynomial(2, 4, 15)", [], 999);
assert_checkerror("gf_minimal_polynomial(7, 1, 7)", [], 999);
