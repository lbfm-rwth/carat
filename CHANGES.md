This file describes changes in various CARAT versions.

* X.Y (YYY-MM-DD)

  - CARAT has been moved to GitHub.
  - Many bugs were fixed, including crashes, wrong results, data corruption
    and purely cosmetic ones
  - The build system was replaced with a completely new one, which made it
    easy to enable additional compiler warnings, which helped track down bugs
  - The C source code for CARAT underwent a major overhaul; constructs
    from old C versions (so called K&R C) were modernized; countless compiler
    warnings were fixed; most turned out harmless, but a seizable number still
    lead us to discover and fix various bugs

* 2.0 (2003-04-11)

  - The new function `Graph` calculates the graph of
    inclusions of a geometric class of space groups.
  - The new function `KSubgroups` calculates maximal klassengleich subgroups.
  - The new function `KSupergroups` calculates minimal klassengleich supergroups
  - The new function `TSubgroups` calculates maximal translationengleich subgroups.
  - The new function `TSupergroups` calculates minimal translationengleich supergroups.
  - Bugfixes have been made for `Z_equiv`, `Name` and `Presentation`.
  - `Extract -r -D` writes space groups directly into files.
  - `Conv` now also converts groups and converts to TeX. 
  - `Name -c` gives the name of a (space-)group in a short form. The
    Hermann-Mauguin-Symbols can be calculated in dimension 2 and 3 now 
    (`Name -M`).
  - `Extensions/Vectorsystems -S` writes representatives of 
    the affine classes in files.

* ...
