Variours Windows-specific installation stuff, used by
../cgal_config.bat

make/			platform-specific makefile and batchfile templates
			in particular, 
make/$(cc)/makefile	is the makefile template obtained from
			preconfigured by Cygwin/install_cgal makefiles
			   and
make/$(cc)/(no)leda/*.bat  contain batchfile templates with lists of
			   demos and examples to compile/link.
		
bin/		binaries used for installation; rebuilt from src/ 
		by cgal_config.bat if necessary 
src/		sources for binaries in bin/
demos/*/makefile.mak	    makefiles for demos
examples/*/makefile.mak	    makefiles for examples

--------------------------------------------------------
Generated (other than by the compiler/linker) files:

make_lib.bat		compiles CGAL.lib
make_demos.bat		compiles demos
make_examples.bat	compiles examples

makefile.mak		main makefile

--------------------------------------------------------
Other files: in

../include/CGAL/config/msvc/
../include/CGAL/config/bcc/

