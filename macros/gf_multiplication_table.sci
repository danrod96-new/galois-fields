// T = gf_multiplication_table(field)
//   Multiplication table of the prime field GF(p): T(i, j) = a * b with
//   a = i - 1 and b = j - 1 (both in 0 .. p - 1).
function T = gf_multiplication_table(field)
    gf_check_field(field, "gf_multiplication_table");
    p = field(1);
    if field(2) <> 1 then
        error(msprintf("%s: Only available for prime fields GF(p).", "gf_multiplication_table"));
    end
    values = 0:p - 1;
    count = size(values, "*");
    T = gf_multiply(field, values' * ones(1, count), ones(count, 1) * values);
endfunction
