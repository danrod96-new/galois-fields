// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// F = gf_primitive_polynomials(p, n)
//   All monic primitive polynomials of degree n over GF(p), one per row
//   (n + 1 coefficients, constant term first).
//   This tests every candidate and is meant for small fields.
function F = gf_primitive_polynomials(p, n)
    fname = "gf_primitive_polynomials";
    // Validates p and n (and that p^n is supported).
    gf_default_primitive_polynomial(p, n);

    if n == 1 then
        field = gf_create_field(p, 1);
        roots = gf_primitive_elements(field);
        F = [gf_mod(-roots(:), p), ones(size(roots, "*"), 1)];
        return
    end

    q = p ^ n;
    F = [];
    for k = 1:q - 1
        coefficients = modulo(floor(k ./ (p .^ (0:n-1))), p);
        if coefficients(1) == 0 then
            continue
        end
        candidate = [coefficients, 1];
        if gf_is_primitive_polynomial(candidate, p, n) then
            F($ + 1, :) = candidate;
        end
    end
endfunction
