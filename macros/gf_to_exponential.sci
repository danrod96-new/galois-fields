// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// k = gf_to_exponential(field, elements [, format])
//   Exponents k such that each element (one per row) equals alpha^k, alpha
//   being a root of the primitive polynomial (extension fields, n > 1).
//   The zero element has exponent -1.
//   format: "e" (default) returns the exponents; "c" returns strings
//   "x^k" ("0" for the zero element).
function k = gf_to_exponential(field, elements, format)
    fname = "gf_to_exponential";
    gf_check_field(field, fname);
    if argn(2) < 3 then
        format = "e";
    end
    p = field(1);
    n = field(2);
    if n < 2 then
        error(msprintf("%s: Only available for extension fields GF(p^n) with n > 1.", fname));
    end
    if type(elements) <> 1 | size(elements, "c") <> n then
        error(msprintf("%s: Wrong size for input argument #%d: %d columns expected (one element per row).", fname, 2, n));
    end
    if or(elements < 0 | elements >= p | elements <> floor(elements)) then
        error(msprintf("%s: Wrong value for input argument #%d: Integers in [0, %d] expected.", fname, 2, p - 1));
    end
    if and(format <> ["e", "c"]) then
        error(msprintf("%s: Wrong value for input argument #%d: ""e"" or ""c"" expected.", fname, 3));
    end

    // Rows of the element table are 0, alpha^0, alpha^1, ...: row r holds
    // alpha^(r - 2). Compare rows through their base-p codes.
    weights = (p .^ (0:n-1))';
    [found, row] = members(elements * weights, field(4) * weights);
    k = row - 2;

    if format == "c" then
        exponents = k;
        k = "x^" + string(exponents);
        k(exponents == -1) = "0";
    end
endfunction
