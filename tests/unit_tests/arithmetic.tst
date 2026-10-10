// Unit tests: arithmetic between field elements.
// <-- NO CHECK REF -->
// <-- CLI SHELL MODE -->

// Prime field GF(7): tables against plain modular arithmetic
f7 = gf_create_field(7, 1);
v = 0:6;
assert_checkequal(gf_addition_table(f7), modulo(v' * ones(1, 7) + ones(7, 1) * v, 7));
assert_checkequal(gf_subtraction_table(f7), pmodulo(v' * ones(1, 7) - ones(7, 1) * v, 7));
assert_checkequal(gf_multiplication_table(f7), modulo(v' * v, 7));
D = gf_division_table(f7);
assert_checkequal(size(D), [6 6]);
assert_checkequal(modulo(D .* (ones(6, 1) * (1:6)), 7), (1:6)' * ones(1, 6));
assert_checkequal(gf_add(f7, 5, 4), 2);
assert_checkequal(gf_subtract(f7, 2, 5), 4);
assert_checkequal(gf_multiply(f7, 3, 5), 1);
assert_checkequal(gf_divide(f7, 1, 3), 5);
assert_checkequal(gf_power(f7, 3, 6), 1);
assert_checkequal(gf_power(f7, 3, 1e12), gf_power(f7, 3, modulo(1e12, 6)));
assert_checkequal(gf_power(f7, 0, 5), 0);
assert_checkequal(gf_multiply(f7, [1 2; 3 4], 2), [2 4; 6 1]);  // one element against many

// GF(2^4): alpha^i * alpha^j = alpha^(i + j)
f = gf_create_field(2, 4);
E = f(4);
nonzero = E(2:16, :);                 // alpha^0 .. alpha^14
for i = 0:14
    expected = E(2 + modulo(i + (0:14), 15), :);
    assert_checkequal(gf_multiply(f, nonzero(i + 1, :), nonzero), expected);
end

// Division undoes multiplication, every element has an inverse
one = [1 0 0 0];
for i = 1:15
    a = nonzero(i, :);
    assert_checkequal(gf_multiply(f, a, gf_divide(f, one, a)), one);
    assert_checkequal(gf_divide(f, gf_multiply(f, a, nonzero), nonzero), ones(15, 1) * a);
end

// Addition is coefficient-wise modulo p; subtraction undoes it
a = [1 0 1 1];
b = [0 1 1 0];
assert_checkequal(gf_add(f, a, b), [1 1 0 1]);
assert_checkequal(gf_subtract(f, gf_add(f, a, b), b), a);

// Powers
alpha = E(3, :);
assert_checkequal(gf_power(f, alpha, 15), one);
assert_checkequal(gf_power(f, alpha, 4), E(6, :));
assert_checkequal(gf_power(f, alpha, 0), one);
assert_checkequal(gf_power(f, alpha, 15 * 1000 + 7), E(9, :));
assert_checkequal(gf_power(f, [0 0 0 0], 15), [0 0 0 0]);   // 0^(q - 1) is 0, not 1
assert_checkequal(gf_power(f, nonzero, 15), ones(15, 1) * one);

// GF(3^2) with a different primitive polynomial (x^2 + 2x + 2)
g = gf_create_field(3, 2, [2 2 1]);
G = g(4);
for i = 0:7
    assert_checkequal(gf_power(g, G(3, :), i), G(2 + i, :));
end
assert_checkequal(gf_multiply(g, G(4, :), G(9, :)), G(2 + modulo(2 + 7, 8), :));

// Errors
assert_checkerror("gf_divide(f, one, [0 0 0 0])", [], 999);
assert_checkerror("gf_divide(f7, 3, 0)", [], 999);
assert_checkerror("gf_power(f, [0 0 0 0], 0)", [], 999);
assert_checkerror("gf_power(f, alpha, -1)", [], 999);
assert_checkerror("gf_add(f, [1 0 1], [1 0 0 0])", [], 999);      // wrong length
assert_checkerror("gf_add(f, [2 0 0 0], [1 0 0 0])", [], 999);    // coefficient >= p
assert_checkerror("gf_add(f7, 7, 1)", [], 999);
assert_checkerror("gf_multiply(f7, [1 2], [1 2 3])", [], 999);
assert_checkerror("gf_add(42, 1, 1)", [], 999);
assert_checkerror("gf_multiplication_table(f)", [], 10000);       // prime fields only
