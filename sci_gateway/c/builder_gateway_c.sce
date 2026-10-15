function builder_gw_c()
    gateway_path = get_absolute_file_path();
    CFLAGS = ilib_include_flag([gateway_path, gateway_path + "../../src/c"]);

    functions = [
        "gf_create_field",                 "sci_gf_create_field",                 "csci6";
        "gf_order",                        "sci_gf_order",                        "csci6";
        "gf_default_primitive_polynomial", "sci_gf_default_primitive_polynomial", "csci6";
        "gf_is_primitive_polynomial",      "sci_gf_is_primitive_polynomial",      "csci6";
        "gf_minimal_polynomial",           "sci_gf_minimal_polynomial",           "csci6";
        "gf_add",                          "sci_gf_add",                          "csci6";
        "gf_subtract",                     "sci_gf_subtract",                     "csci6";
        "gf_multiply",                     "sci_gf_multiply",                     "csci6";
        "gf_divide",                       "sci_gf_divide",                       "csci6";
        "gf_power",                        "sci_gf_power",                        "csci6";
        "gf_poly_add",                     "sci_gf_poly_add",                     "csci6";
        "gf_poly_subtract",                "sci_gf_poly_subtract",                "csci6";
        "gf_poly_multiply",                "sci_gf_poly_multiply",                "csci6";
        "gf_poly_divide",                  "sci_gf_poly_divide",                  "csci6";
        "gf_poly_power",                   "sci_gf_poly_power",                   "csci6";
        "gf_is_prime",                     "sci_gf_is_prime",                     "csci6";
        "gf_mod",                          "sci_gf_mod",                          "csci6"
    ];

    files = ["gw_helpers.c", "sci_gf_field.c", "sci_gf_elements.c", "sci_gf_polynomials.c"];

    tbx_build_gateway("gateway_galois_fields", functions, files, gateway_path, ..
                      "../../src/c/libgalois_fields", "", CFLAGS);
endfunction

builder_gw_c();
clear builder_gw_c;
