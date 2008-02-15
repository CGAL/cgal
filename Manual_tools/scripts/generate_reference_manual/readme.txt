generate_reference_manual
*************************

generate_reference_manual automatically generates (most of) the Reference Manual of a CGAL package. It runs Doxygen to generate a Latex documentation from C++ source code, then converts it to CGAL Manual format.
If the C++ source code is commented, this tool generates the whole Reference Manual of the package (except the introduction page). If not, it generates a raw latex documentation that you have to complete manually.

Note that this tool adds or updates automatic documentation and does *not* remove the documentation manually written. Tags %START-AUTO and %END-AUTO clearly deliminate the automatic documentation.

Installation:

This is a shell script, thus for Unix only.

generate_reference_manual and copy_doxygen_latex_doc must be in the PATH.
generate_reference_manual_Doxyfile must be in the same folder as generate_reference_manual.

generate_reference_manual requires Perl (any recent version) and Doxygen >= 1.4.

Usage:

generate_reference_manual [options] /path/to/package/root
    -h, --help      Print this help
    -d, --debug     Turn on debug traces

Typical scenario:

1) Create the package's Reference Manual as described in CGAL Developer Manual (using cc_ref_wizard).
2) Run generate_reference_manual.
3) Run cgal_manual.
4) Look at the generated documentation in the Reference Manual .tex files and in the PS/PDF/HTML final Reference manual. In order to complete it, you may:
    a) Comment the C++ source code using Doxygen conventions and run generate_reference_manual again.
    b) Write extra documentation in .tex files outside of automatic sections.
5) Goto point 3) until the Reference Manual is complete.

Tips:

* To generate a concept's documentation, you may create a fake C++ class named after the concept anywhere in the package's folders tree (use the dont_submit file to ignore it in CGAL releases).
* generate_reference_manual generates a documentation for public and protected items only.
* "Words surrounded by double quotes" are emphasized. Words containing an underscore are automatically emphasized.
* A comment starting by 'Concept:' is converted to a \ccIsModel paragraph.
* A comment starting by 'Sub-concept:' is converted to a \ccRefines paragraph.
* A comment starting by 'Models:' is converted to a \ccHasModels paragraph.
* A comment starting by 'Design pattern:' is converted to a \ccHeading{Design Pattern} paragraph.
* A comment starting by 'Template parameters:' is converted to \ccParameters paragraph.
* 2 comments '@cond SKIP_IN_MANUAL' and '@endcond' delimit a section ignored by generate_reference_manual.

Known bugs:

* So far, this tool documents only concepts, classes, structs and functions.

Contact:

Please contact Laurent Saboret <Laurent.Saboret@sophia.inria.fr> if you discover a bug or wish to request for an enhancement.
