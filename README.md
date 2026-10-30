# galois_fields — finite field arithmetic for Scilab

`galois_fields` is a Scilab toolbox for working with finite (Galois) fields
GF(p) and GF(p^n). It covers:

- building fields from a primitive polynomial;
- element arithmetic: add, subtract, multiply, divide and power;
- arithmetic on polynomials over GF(p);
- primitive and minimal polynomials;
- primitive elements;
- operation tables;
- conversion between polynomial and exponential (α^k) notation.

## Requirements

- Scilab 6.1 or later. Tested with Scilab 2024.0.0. It uses the Scilab 6 C
  gateway API (`api_scilab.h`, "csci6" gateways), which is the current API.
- A C compiler that Scilab can use (`haveacompiler()` must return `T`): gcc
  on Linux, Visual Studio Build Tools on Windows, or clang on macOS. See
  [Compiling on each system](#compiling-on-each-system) for the setup on each one.

## Building and loading

From the toolbox folder:

```scilab
exec builder.sce   // compiles the C library and gateways, builds macros and help pages
exec loader.sce    // loads the toolbox into the current session
```

The help pages are compiled into `jar/` with Scilab's Java tools, so run
`builder.sce` from the normal Scilab desktop or `scilab -nw`. With
`scilab-cli` everything else is built, and the help is skipped with a message.
Once the toolbox is loaded, `help gf_create_field` (or any other function name)
opens its page, and the "Galois Fields" book appears in the help browser.

`exec unloader.sce` unloads it, and `exec cleaner.sce` removes the build
products. To load it automatically, add `exec("<path>/loader.sce", -1);` to your
startup file (`SCIHOME/.scilab` or `SCIHOME/scilab.ini`).

Run the test suite with:

```scilab
test_run(pwd(), [], ["no_check_error_output"])
```

The demos ("Arithmetic in GF(2^4)" and "Hamming (7,4) error control") appear
in the demo browser (`demo_gui()`) once the toolbox is loaded.

## Compiling on each system

Compiled binaries for Linux, Windows and macOS are available from ATOMS
(`atomsInstall("galois_fields")`), which needs no compiler. Follow these steps
only to build the toolbox from source.

The steps are the same everywhere once a compiler is installed:

1. Install Scilab 6.1 or later from <https://www.scilab.org/download>. Use the
   full desktop version, since the help is built with its Java tools.
2. Install a C/C++ compiler for your system (see below). Then restart Scilab
   and check that it finds the compiler:

   ```scilab
   haveacompiler()    // must return T
   ```

3. Build, load and test the toolbox from its folder:

   ```scilab
   cd("<path>/galois_fields")
   exec("builder.sce", -1)
   exec("loader.sce", -1)
   test_run(pwd(), [], ["no_check_error_output"])
   ```

4. Optionally, pack an ATOMS binary archive for this system with
   `exec("make_atoms_archive.sce", -1)`. It checks that every built file is
   present and writes `galois_fields-<version>-<os>.tar.gz` (or `.zip` on
   Windows) next to the toolbox folder.

Binaries only work on the system and Scilab version they were built with, so
repeat these steps on each system and for each Scilab version you target.

### Linux

1. Install gcc, g++ and make:

   ```sh
   sudo apt install build-essential        # Debian, Ubuntu, Mint
   sudo dnf install gcc gcc-c++ make       # Fedora, RHEL, Rocky
   sudo pacman -S base-devel               # Arch, Manjaro
   ```

2. If you use your distribution's Scilab package rather than the scilab.org
   download, also install its development headers when they are packaged
   separately (for example `scilab-include` on Debian and Ubuntu).
3. Follow the common steps. The build creates `src/c/libgalois_fields.so` and
   `sci_gateway/c/libgateway_galois_fields.so`.

A library built on a recent distribution may not load on older ones, because
it needs at least the glibc version it was built against. To publish a binary
that works on most systems, build it on an older, long-term-support release.

### Windows (64-bit)

1. Install the free **Build Tools for Visual Studio** from
   <https://visualstudio.microsoft.com/downloads/> (under "Tools for Visual
   Studio"). In the installer, select the **Desktop development with C++**
   workload, which includes the MSVC compiler and the Windows SDK.
   - Scilab 2023 and later: Visual Studio 2022.
   - Scilab 6.1: Visual Studio 2019, because 6.1 does not detect newer
     versions.
2. Restart Scilab and check `haveacompiler()`. If it returns `F`, try the
   other Visual Studio version, or on Scilab 6.x install MinGW from ATOMS with
   `atomsInstall("mingw")` and restart Scilab.
3. Keep the toolbox in a folder whose path has no spaces or accented letters,
   such as `C:\scilab\galois_fields`, since some compilers fail on such paths.
4. Follow the common steps. The build creates `src/c/libgalois_fields.dll` and
   `sci_gateway/c/libgateway_galois_fields.dll`.

### macOS

1. Install Apple's command line tools, which provide the clang compiler:

   ```sh
   xcode-select --install
   ```

   The full Xcode app is not needed.
2. Install the Scilab build that matches your Mac: Apple Silicon (arm64) or
   Intel (x86_64). The first time, open Scilab with right-click → Open if
   macOS blocks it.
3. Follow the common steps. The build creates `src/c/libgalois_fields.dylib`
   and `sci_gateway/c/libgateway_galois_fields.dylib`.

The libraries match the architecture of the Scilab that built them. To publish
for both kinds of Mac, build once with each Scilab build.

## Data conventions

- **Polynomials over GF(p)** are row vectors of coefficients in `0 .. p-1`,
  constant term first: `[1 1 0 0 1]` is 1 + x + x⁴. Results are trimmed to
  their degree, and the zero polynomial is `0`.
- **Fields** are the lists returned by `gf_create_field`:
  `list(p, n, primitive_polynomial, elements)`.
  - For n = 1, `elements` is `0:p-1`.
  - For n > 1, `elements` is a p^n × n matrix in exponential order: row 1 is
    0, and row k + 2 is α^k, where α is a root of the primitive polynomial.
- **Elements** of GF(p) are integers. Elements of GF(p^n), n > 1, are rows of
  n coefficients (polynomials in α of degree < n, constant term first).
  Element functions accept several elements at once (a matrix for n = 1, or
  one element per row for n > 1). A single element is combined with each of
  the others.
- **Limits:** p^n ≤ 2^20 (1 048 576) and n ≤ 20.

## Quick example

```scilab
F = gf_create_field(2, 4);           // GF(16); default primitive poly 1 + x + x^4
F(3)                                 // [1 1 0 0 1]
gf_order(F)                          // 16
a = [0 1 0 0];                       // α
b = [1 1 0 0];                       // 1 + α = α^4
c = gf_multiply(F, a, b)             // [0 1 1 0] = α + α^2
gf_to_exponential(F, c)              // 5, i.e. α^5
gf_power(F, a, 15)                   // [1 0 0 0] = 1
gf_minimal_polynomial(F, 3)          // [1 1 1 1 1] = 1 + x + x^2 + x^3 + x^4
[r, q] = gf_poly_divide([1 0 1 1], [1 1], 2)   // r = 1, q = [0 0 1]
```

## Function reference

### Fields and polynomials (C gateways)

| Function | Description |
|---|---|
| `F = gf_create_field(p, n [, f])` | Field GF(p^n). `f` (n+1 coefficients) must be primitive; if it is omitted, the default primitive polynomial is used. |
| `q = gf_order(F)` | Number of elements, p^n. |
| `f = gf_default_primitive_polynomial(p, n)` | Default primitive polynomial. For n = 1 this is x − g, with g the smallest primitive root; for n > 1 it is the first monic primitive polynomial. |
| `tf = gf_is_primitive_polynomial(f, p, n)` | `T` if `f` (degree n) is primitive over GF(p). |
| `m = gf_minimal_polynomial(p, n, e)` or `gf_minimal_polynomial(F, e)` | Minimal polynomial over GF(p). For n = 1, `e` is the element; for n > 1, it is the exponent k of α^k (−1 for zero). |
| `c = gf_add(F, a, b)`, `gf_subtract`, `gf_multiply`, `gf_divide` | Element arithmetic. Dividing by zero is an error. |
| `c = gf_power(F, a, e)` | a^e for an integer e ≥ 0 (0^0 is an error). |
| `c = gf_poly_add(a, b, p)`, `gf_poly_subtract`, `gf_poly_multiply` | Polynomial arithmetic over GF(p). The arguments may have different lengths. |
| `[r, q] = gf_poly_divide(a, b, p)` | Division with remainder: a = q·b + r. |
| `c = gf_poly_power(a, e, p)` | a^e (result degree ≤ 100 000). |
| `tf = gf_is_prime(x)` | Element-wise primality test. |
| `r = gf_mod(x, m)` | Non-negative remainder of integers modulo m. |

### Macros

| Function | Description |
|---|---|
| `p = gf_characteristic(F)`, `n = gf_degree(F)` | Field parameters. |
| `f = gf_primitive_polynomial(F)` | The field's primitive polynomial as a Scilab `poly`. |
| `e = gf_elements(F)` | All elements as a column of Scilab polynomials. |
| `P = gf_to_polynomial(coefficients)` | Coefficient rows → Scilab polynomials. |
| `k = gf_to_exponential(F, elements [, "e" \| "c"])` | For n > 1: exponent k of α^k for each row (zero gives −1), or strings such as `"x^5"` / `"0"` with `"c"`. |
| `g = gf_primitive_elements(F)` | All primitive elements. |
| `F = gf_primitive_polynomials(p, n)` | All monic primitive polynomials of degree n, one per row. |
| `gf_addition_table(F)`, `gf_subtraction_table`, `gf_multiplication_table`, `gf_division_table` | Operation tables for prime fields (n = 1). |
| `gf_is_zero(a)`, `gf_is_one(a)`, `gf_are_equal(a, b)` | Small helpers on coefficient arrays. |

Every function has a help page with its arguments, a description and a
working example; start with `help galois_fields` for an overview. The page
sources are in `help/en_US` and in `help/es_ES` (Spanish).

## Summary of changes

- Invalid arguments raise real Scilab errors, which you can catch with
  `try`/`catch`. The old functions displayed a message and returned 0.
- The functions no longer print informational text.
- `gf_create_field(p, n)` picks a default primitive polynomial for every n,
  not only n = 1.
- `gf_is_primitive_polynomial` and `gf_is_prime` return booleans.
- Element operations accept several elements at once.
- Polynomial addition and subtraction accept operands of different lengths
  and trim the result.
- `gf_primitive_elements` no longer has the p ≤ 107 limit.
- Fixed: raising a non-zero element to a multiple of q − 1 now gives 1.
- The Hamming demo was rewritten. It no longer asks for input, and it shows
  how each single-bit error is located from its syndrome.

## Not included yet

- **Translations of the help.** Only English and Spanish pages exist. Pages in another
  language go in a sibling folder such as `help/fr_FR`, which `builder.sce`
  picks up automatically.

## License

This toolbox is released under the terms of the CeCILL license :
http://www.cecill.info/index.en.html

## Author

Copyright (C) 2009-2026  Daniel E Rodriguez — Universidad Industrial de Santander (UIS).
