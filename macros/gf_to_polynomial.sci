// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// P = gf_to_polynomial(coefficients)
//   Converts coefficient rows (constant term first) to Scilab polynomials
//   in x, one per row.
function P = gf_to_polynomial(coefficients)
    if type(coefficients) <> 1 | ~isreal(coefficients) then
        error(msprintf("%s: Wrong type for input argument #%d: A real matrix expected.", "gf_to_polynomial", 1));
    end
    rows = size(coefficients, "r");
    P = poly(0, "x", "coeff");
    for i = 1:rows
        P(i, 1) = poly(coefficients(i, :), "x", "coeff");
    end
endfunction
