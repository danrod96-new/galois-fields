// Packs the built galois_fields toolbox into an ATOMS binary archive for the
// current OS. Run it from Scilab after exec("builder.sce", -1):
//     exec("make_atoms_archive.sce", -1)
// Result: galois_fields-<version>-<os>.zip (Windows) or .tar.gz (Linux, macOS),
// written next to the toolbox folder.
mode(-1);
lines(0);

function make_atoms_archive()
    TOOLBOX_NAME = "galois_fields";
    root = get_absolute_file_path();

    // Version from DESCRIPTION
    desc = mgetl(root + "DESCRIPTION");
    v = desc(grep(desc, "/^Version:/", "r"));
    version = stripblanks(part(v(1), length("Version:") + 1:length(v(1))));

    // Platform
    os = getos();
    ext = getdynlibext();
    select os
    case "Windows" then osname = "windows64";
    case "Darwin"  then osname = "macos";
    else                osname = "linux64";
    end

    // Files that must be in the archive; the build must have created them.
    required = [
        "DESCRIPTION"
        "license.txt"
        "changelog.txt"
        "README.md"
        "loader.sce"
        "unloader.sce"
        "etc/galois_fields.start"
        "etc/galois_fields.quit"
        "macros/lib"
        "src/c/loader.sce"
        "src/c/libgalois_fields" + ext
        "sci_gateway/loader_gateway.sce"
        "sci_gateway/c/loader.sce"
        "sci_gateway/c/libgateway_galois_fields" + ext
        "jar/scilab_en_US_help.jar"
        "jar/scilab_es_ES_help.jar"
        "demos/galois_fields.dem.gateway.sce"
    ];
    missing = [];
    for f = required'
        if ~isfile(root + f) then
            missing = [missing; f];
        end
    end
    if missing <> [] then
        error(msprintf("%s: missing, run builder.sce from the Scilab desktop first:\n  %s\n", ..
            TOOLBOX_NAME, strcat(missing, "\n  ")));
    end

    // Stage a clean copy in TMPDIR/galois_fields
    stage_root = TMPDIR + filesep() + "atoms_stage";
    stage = stage_root + filesep() + TOOLBOX_NAME;
    if isdir(stage_root) then
        rmdir(stage_root, "s");
    end
    mkdir(stage_root);
    mkdir(stage);
    for d = ["etc" "macros" "src" "src/c" "sci_gateway" "sci_gateway/c" "jar" "demos"]
        mkdir(stage + filesep() + d);
    end

    for f = required'
        copyfile(root + f, stage + filesep() + f);
    end
    // Compiled and source macros (sources help users read the code)
    for pat = ["*.bin" "*.sci"]
        for f = listfiles(root + "macros/" + pat)'
            copyfile(f, stage + "/macros/");
        end
    end
    // All demos
    for f = listfiles(root + "demos/*.sce")'
        copyfile(f, stage + "/demos/");
    end
    // Unit tests (sources only, not .dia outputs)
    if isdir(root + "tests/unit_tests") then
        mkdir(stage + "/tests");
        mkdir(stage + "/tests/unit_tests");
        for f = listfiles(root + "tests/unit_tests/*.tst")'
            copyfile(f, stage + "/tests/unit_tests/");
        end
    end

    // Compress next to the toolbox folder
    parent = fullpath(root + "..");
    base = TOOLBOX_NAME + "-" + version + "-" + osname;
    if os == "Windows" then
        archive = parent + filesep() + base + ".zip";
        if isfile(archive) then deletefile(archive); end
        cmd = "powershell -NoProfile -Command ""Compress-Archive -Path ''" + stage + ..
              "'' -DestinationPath ''" + archive + "''""";
    else
        archive = parent + filesep() + base + ".tar.gz";
        cmd = "tar -czf """ + archive + """ -C """ + stage_root + """ " + TOOLBOX_NAME;
    end
    if host(cmd) <> 0 then
        error(msprintf("%s: compression failed: %s\n", TOOLBOX_NAME, cmd));
    end
    rmdir(stage_root, "s");
    mprintf("%s: created %s\n", TOOLBOX_NAME, archive);
endfunction

make_atoms_archive();
clear make_atoms_archive;
