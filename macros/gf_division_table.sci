// T = gf_division_table(field)
//   Division table of the prime field GF(p): T(i, j) = a / b with
//   a = i and b = j (both in 1 .. p - 1).
function T = gf_division_table(field)
    gf_check_field(field, "gf_division_table");
    p = field(1);
    if field(2) <> 1 then
        error(msprintf("%s: Only available for prime fields GF(p).", "gf_division_table"));
    end
    values = 1:p - 1;
    count = size(values, "*");
    T = gf_divide(field, values' * ones(1, count), ones(count, 1) * values);
endfunction
