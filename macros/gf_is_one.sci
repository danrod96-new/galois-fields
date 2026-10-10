// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// tf = gf_is_one(a)
//   %t if a is the unit element: 1 followed only by zeros.
function tf = gf_is_one(a)
    a = a(:);
    tf = a(1) == 1 & and(a(2:$) == 0);
endfunction
