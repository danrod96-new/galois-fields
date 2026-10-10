// Unit tests: representations, primitive elements and helpers.
// <-- NO CHECK REF -->
// <-- CLI SHELL MODE -->

f = gf_create_field(2, 4);
E = f(4);

// Exponential representation: row r of the table is alpha^(r - 2)
assert_checkequal(gf_to_exponential(f, E), (-1:14)');
assert_checkequal(gf_to_exponential(f, [1 1 0 0]), 4);
s = gf_to_exponential(f, E(1:3, :), "c");
assert_checkequal(s, ["0"; "x^0"; "x^1"]);
assert_checkerror("gf_to_exponential(gf_create_field(7, 1), 3)", [], 10000);
assert_checkerror("gf_to_exponential(f, [2 0 0 0])", [], 10000);
assert_checkerror("gf_to_exponential(f, E, ""z"")", [], 10000);

// Polynomial representation
x = poly(0, "x");
assert_checkequal(gf_primitive_polynomial(f), 1 + x + x^4);
assert_checkequal(gf_to_polynomial([1 1 0 1]), 1 + x + x^3);
P = gf_elements(f);
assert_checkequal(size(P), [16 1]);
assert_checkequal(P(6), 1 + x);
P7 = gf_elements(gf_create_field(7, 1));
assert_checkequal(size(P7), [7 1]);
assert_checkequal(coeff(P7(4)), 3);

// Primitive elements
assert_checkequal(gf_primitive_elements(gf_create_field(7, 1)), [3 5]);
assert_checkequal(gf_primitive_elements(gf_create_field(2, 1)), 1);
assert_checkequal(gf_primitive_elements(gf_create_field(23, 1)), [5 7 10 11 14 15 17 19 20 21]);
assert_checkequal(size(gf_primitive_elements(gf_create_field(113, 1)), "*"), 48);  // phi(112), beyond the old table
prim = gf_primitive_elements(f);
assert_checkequal(size(prim, "r"), 8);                                              // phi(15)
for i = 1:8
    // a primitive element has order 15: its powers 15/5 and 15/3 are not 1
    for k = [3 5]
        assert_checkfalse(gf_is_one(gf_power(f, prim(i, :), k)));
    end
end

// Helpers
assert_checktrue(gf_is_zero([0 0 0]));
assert_checkfalse(gf_is_zero([0 1 0]));
assert_checktrue(gf_is_one([1 0 0]));
assert_checkfalse(gf_is_one([0 1 0]));
assert_checktrue(gf_are_equal([1 2], [1 2]));
assert_checkfalse(gf_are_equal([1 2], [1 2 0]));
