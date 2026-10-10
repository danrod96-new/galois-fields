// Demos of the galois_fields toolbox.
function subdemolist = demo_gateway()
    demopath = get_absolute_file_path();
    subdemolist = ["Arithmetic in GF(2^4)",       "field_arithmetic.dem.sce";
                   "Hamming (7,4) error control", "hamming_code.dem.sce"];
    subdemolist(:, 2) = demopath + subdemolist(:, 2);
endfunction

subdemolist = demo_gateway();
clear demo_gateway;
