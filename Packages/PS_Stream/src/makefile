
default: io_lib

      

# Please adapt it to your environment

# (1) The Compiler  
# 
#     Pay attention not to put a blank  behind the compiler name!!!

# for GNU Compiler

#      CC   = g++

# for SGI and SUN Pro Compiler
CC = CC

# (2) The operating system of the machine where you are compiling.
#
#     Issue the command "sh system" in the CGAL home directory
#     to find out what OS you are using. 

#     It is a good idea to ut all these lines in comments and to
#     define an environment variable in your .profile or .login
#
#     Pay attention not to put a blank behind the OS name!!!
#
#CGAL_OSTYPE =irix5
#CGAL_OSTYPE =sun4OS4
CGAL_OSTYPE =sun4OS5



# (3) The Home Directory of LEDA
#
#     It is a good idea to create a symbolic link to the distribution
#
LEDA_DIR = $(HOME)/Leda



# (4) The Home Directory of CGAL
#
#     It is a good idea to create a symbolic link to the distribution
#
CGAL_DIR = $(HOME)/Cgal


# If the makefiles for Cgal are not under $(CGAL_DIR) you
# have to change the following:

CGAL_MAKE_DIR = $(CGAL_DIR)/MAKE

# If the header files and the libraries are not installed in 
# subdirectories of the LEDA/CGAL root files, you may have to
# change the following too:
 

CGAL_INCL_DIR = $(CGAL_DIR)/include
CGAL_LIB_DIR = $(CGAL_DIR)/lib/$(CGAL_OSTYPE)_$(CC)


LEDA_INCL_DIR = $(LEDA_DIR)/incl
LEDA_LIB_DIR = $(LEDA_DIR)/lib/$(CGAL_OSTYPE)_$(CC)




#..........................................................................
# Don't edit the following line
include $(CGAL_MAKE_DIR)/$(CGAL_OSTYPE)_$(CC)/application
#..........................................................................


      
      
CFLAGS = -g  $(CGAL_CFLAGS)
LDFLAGS = $(CGAL_LDFLAGS)

all: PS_Stream

PS_Stream: Postscript_stream.o postscript.o
	$(CC) -o PS_Stream Postscript_stream.o postscript.o $(LDFLAGS)
Demo: Postscript_stream.o Demo.o
	$(CC) -o Demo Postscript_stream.o Demo.o $(LDFLAGS)
.C.o:
	@echo "compiling "  $(notdir $<)
	$(CC)  $(CFLAGS)  -c $<

clean:
	rm *.o

ptclean:
	ptclean
