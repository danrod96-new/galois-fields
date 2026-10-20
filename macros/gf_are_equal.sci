// tf = gf_are_equal(a, b)
//   T if a and b have the same size and the same coefficients.
function tf = gf_are_equal(a, b)
    tf = and(size(a) == size(b));
    if tf then
        tf = and(a == b);
    end
endfunction
