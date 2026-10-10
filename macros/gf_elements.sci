// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// e = gf_elements(field)
//   Elements of the field as a column of Scilab polynomials in x: 0 first,
//   then alpha^0, alpha^1, ..., alpha^(p^n - 2) for n > 1, or 0 .. p - 1
//   for n = 1. field(4) holds the same elements as coefficient rows.
function e = gf_elements(field)
    gf_check_field(field, "gf_elements");
    table = field(4);
    if field(2) == 1 then
        table = table(:);
    end
    e = gf_to_polynomial(table);
endfunction
