There are two ways to compile: inside the CGAL tree or outside it. Clearly,
it is easier to do it the second way, but compiling inside CGAL will give you
all the advantages of the install script.

To compile inside the tree, it is necessary to copy all the files to the CGAL
root directory, and install it in the normal way, using the provided modified
script install_cgal.

To compile outside the tree, just go to the directory src/Gbrs and type
"make outside". It will install the library in the CGAL lib directory.
Note that there must be present the environment variables CGAL_MAKEFILE,
RS_INCL, RS_LIB, MPFI_INCL and MPFI_LIB.

Have fun,

--
Luis

