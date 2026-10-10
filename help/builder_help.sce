// Builds the help of the galois_fields toolbox (one jar per language).
here = get_absolute_file_path();
d = dir(here);
d = d.name(d.isdir);
la_LA = d(grep(d, "/[a-z]{2}_[A-Z]{2}/", "r"));
tbx_builder_help_lang(la_LA, here);
clear here d la_LA
