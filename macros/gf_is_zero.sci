// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// tf = gf_is_zero(a)
//   %t if every coefficient of a is 0.
function tf = gf_is_zero(a)
    tf = and(a == 0);
endfunction
