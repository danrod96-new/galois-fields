// Builds the galois_fields toolbox: macros, C library and gateways.
// Run it from Scilab:  exec("builder.sce", -1)
mode(-1);
lines(0);

function main_builder()
    TOOLBOX_NAME  = "galois_fields";
    TOOLBOX_TITLE = "Galois Fields";
    toolbox_dir   = get_absolute_file_path();

    if ~haveacompiler() then
        error(msprintf("%s: A C compiler is required to build this toolbox.", TOOLBOX_NAME));
    end

    tbx_builder_macros(toolbox_dir);
    tbx_builder_src(toolbox_dir);
    tbx_builder_gateway(toolbox_dir);

    // The help is compiled to jar/ by Java tools, which scilab-cli lacks.
    if getscilabmode() <> "NWNI" then
        tbx_builder_help(toolbox_dir);
    else
        mprintf("%s: Help not built (it needs Scilab with Java, not scilab-cli).\n", TOOLBOX_NAME);
    end

    tbx_build_loader(toolbox_dir);
    tbx_build_cleaner(toolbox_dir);
endfunction

main_builder();
clear main_builder;
