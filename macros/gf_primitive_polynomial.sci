// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// f = gf_primitive_polynomial(field)
//   Primitive polynomial of the field as a Scilab polynomial in x.
function f = gf_primitive_polynomial(field)
    gf_check_field(field, "gf_primitive_polynomial");
    f = poly(field(3), "x", "coeff");
endfunction
