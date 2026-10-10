function buildmacros()
    macros_path = get_absolute_file_path();
    tbx_build_macros("galois_fields", macros_path);
endfunction

buildmacros();
clear buildmacros;
