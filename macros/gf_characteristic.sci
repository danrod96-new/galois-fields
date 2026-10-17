// p = gf_characteristic(field)
//   Characteristic (prime p) of the field GF(p^n).
function p = gf_characteristic(field)
    gf_check_field(field, "gf_characteristic");
    p = field(1);
endfunction
