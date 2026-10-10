// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// n = gf_degree(field)
//   Degree n of the field GF(p^n) over its prime subfield GF(p).
function n = gf_degree(field)
    gf_check_field(field, "gf_degree");
    n = field(2);
endfunction
