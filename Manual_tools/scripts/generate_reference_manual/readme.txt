generate_reference_manual
*************************

generate_reference_manual automatically generates (most of) the Reference Manual of a CGAL package. It runs Doxygen to generate a Latex documentation from C++ source code, then converts it to CGAL Manual format.
If the C++ source code is commented, this tool generates the whole Reference Manual of the package (except the introduction page). If not, it generates a raw latex documentation that you will have to complete manually.

Note that this tool adds or updates automatic documentation and does *not* remove the documentation manually written. Tags %START-AUTO and %END-AUTO clearly deliminate the automatic documentation.

INSTALLATION

This is a shell script, thus for Unix only.

generate_reference_manual and copy_doxygen_latex_doc must be in the PATH.
generate_reference_manual_Doxyfile must be in the same folder as generate_reference_manual.

generate_reference_manual requires Perl (any recent version) and Doxygen (Doxygen 1.4 or 1.5).

USAGE

generate_reference_manual [options] /path/to/package/root
    -h, --help      Print this help
    -f, --force     Insert missing %START-AUTO..%END-AUTO sections

TYPICAL SCENARIO

1) Create the package's Reference Manual as described in CGAL Developer Manual.
2) Run generate_reference_manual --force to create the %START-AUTO..%END-AUTO sections.
3) Run cgal_manual.
4) Look at the generated documentation in the Reference Manual .tex files and in the PS/PDF/HTML final Reference manual. In order to complete it, you may:
    a) Comment the C++ source code using Doxygen conventions.
    b) Write extra documentation in .tex files outside of automatic sections.
    c) Modify automatic sections (make sure to delete the %START-AUTO and %END-AUTO tags).
5) Run generate_reference_manual. Do not use the --force option anymore if you modified automatic sections.
6) Goto point 3) until the Reference Manual is complete.

TIPS

* To generate a concept's documentation, you may create a fake C++ class named after the concept anywhere in the package's folders tree (use the dont_submit file to ignore it in CGAL releases).
* A comment starting by 'Concept:' will be converted to a \ccIsModel paragraph.
* A comment starting by 'Sub-concept:' will be converted to a \ccRefines paragraph.
* A comment starting by 'Models:' will be converted to a \ccHasModels paragraph.
* A comment starting by 'Design pattern:' will be converted to a \ccHeading{Design Pattern} paragraph.
* A comment starting by 'Template parameters:' will be converted to \ccParameters paragraph.
* "Words surrounded by double quotes" will be emphasized.

KNOWNS BUGS

* This tool has been tested on Linux with Doxygen 1.4 and 1.5. It is likely that it will *not* work with other versions of Doxygen.
* So far, this tool documents only concepts, classes, structs and functions.

CONTACT

Please contact Laurent Saboret <Laurent.Saboret@sophia.inria.fr> if you discover a bug or wish to request for an enhancement.
