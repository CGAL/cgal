# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = 
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
#

CC = $(CGAL_CXX)

CXXFLAGS = \
           -I../include \
           $(CGAL_CXXFLAGS) \
           -g 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          -L../lib -lVanill \
          $(CGAL_WINDOW_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

TARGETS = vanilla_test1 vanilla_test2 vanilla_test3

TEST1_OBJECTS = vanilla_test1.o  

TEST2_OBJECTS = vanilla_test2.o 

TEST3_OBJECTS = vanilla_test3.o 

all:   $(TARGETS)

vanilla_test1 : $(TEST1_OBJECTS)
	$(CC) -o vanilla_test1 $(TEST1_OBJECTS) $(LDFLAGS) 

vanilla_test2 : $(TEST2_OBJECTS)
	$(CC) -o vanilla_test2 $(TEST2_OBJECTS) $(LDFLAGS) 

vanilla_test3 : $(TEST3_OBJECTS)
	$(CC) -o vanilla_test3 $(TEST3_OBJECTS) $(LDFLAGS) 

vanilla_test1.o: vanilla_test1.C\
	../include/CEP/Vanilla/Flavored_object.h 

vanilla_test2.o: vanilla_test2.C\
	../include/CEP/Vanilla/Flavored_object.h 

vanilla_test3.o: vanilla_test3.C\
	../include/CEP/Vanilla/Flavored_object.h 

flavor.o: \
	../include/CEP/Vanilla/flavor.h 

clean:
	/bin/rm -f *.o \
        core

really_clean: clean
	/bin/rm -f *.out \
        $(TARGETS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CC) $(CXXFLAGS) -c  $<

