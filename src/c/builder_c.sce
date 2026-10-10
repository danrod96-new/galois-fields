// Builds libgalois_fields, the C library that does the field arithmetic.
function builder_c()
    src_c_path = get_absolute_file_path();
    CFLAGS = ilib_include_flag(src_c_path);
    tbx_build_src("galois_fields", "galois_fields.c", "c", src_c_path, "", "", CFLAGS);
endfunction

builder_c();
clear builder_c;
