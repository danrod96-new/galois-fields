// g = gf_primitive_elements(field)
//   Primitive elements (generators of the multiplicative group) of the field.
//   For n = 1: a row of integers in increasing order.
//   For n > 1: one element per row, alpha^k for every k coprime to p^n - 1.
function g = gf_primitive_elements(field)
    gf_check_field(field, "gf_primitive_elements");
    p = field(1);
    n = field(2);
    q = p ^ n;

    if n == 1 then
        if p == 2 then
            g = 1;
            return
        end
        // The primitive polynomial is c1*x + c0, whose root is -c0/c1.
        f = field(3);
        root = gf_divide(field, gf_mod(-f(1), p), f(2));
        g = [];
        for k = 1:p - 1
            if gcd(int32([k, p - 1])) == 1 then
                g($ + 1) = gf_power(field, root, k);
            end
        end
        g = gsort(g', "g", "i");
    else
        table = field(4);
        g = [];
        for k = 1:q - 2
            if gcd(int32([k, q - 1])) == 1 then
                g($ + 1, :) = table(k + 2, :);
            end
        end
    end
endfunction
