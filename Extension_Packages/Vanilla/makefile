# This is the makefile for compiling and installing the library libVanill.a 

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = 
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                directory where CGAL is installed
#---------------------------------------------------------------------#
#CGAL_HOME=


#---------------------------------------------------------------------#
#                    installation program
#---------------------------------------------------------------------#
INSTALL = cp
INSTALL_FLAGS = -R

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#


library:
	cd src; make

install:
	$(INSTALL) $(INSTALL_FLAGS) include/* $(CGAL_HOME)/include   
	$(INSTALL) $(INSTALL_FLAGS) lib/* $(CGAL_HOME)/lib/$(CGAL_OS_COMPILER) 

clean:
	cd src; make clean
	cd test_suite; make clean

really_clean:
	cd src; make really_clean
	cd test_suite; make really_clean
	rm -f lib/libVanill.a

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CC) $(CXXFLAGS) -c  $<

