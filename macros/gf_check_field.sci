// Copyright (C) Daniel E Rodriguez - UIS.
// Part of the galois_fields toolbox.

// gf_check_field(field, caller)
//   Raises an error, attributed to caller, unless field is a field created
//   by gf_create_field.
function gf_check_field(field, caller)
    ok = typeof(field) == "list" & size(field) == 4;
    if ok then
        try
            gf_order(field);
        catch
            ok = %f;
        end
    end
    if ~ok then
        error(msprintf("%s: Wrong type for input argument #%d: A field created by gf_create_field expected.", caller, 1));
    end
endfunction
