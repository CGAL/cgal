#  Makefile
#  ===================================
#     V 1.00   17.08.1995   Lutz Kettner
#              Project: CGAL merger tool for the specification task
#     V 1.1    22.11.1996   Lutz Kettner
#              Project: tools around the CGAL manual writing task.
#              extract, check, HTML conversion.
#     V 1.2    03.12.1996   Lutz Kettner
#              Tools renamed from CGAL to CC. Relates to cc_manual.sty R2.1.
#
#     $Release:$, $Date$

# =======================================================
# Installation: chose the appropriate compiler and paths:
# =======================================================

# Installation directory (where to put the executables)
# Only necessary if the make target 'install' is used for installation.
#INSTALL_DIR = /usr/bin
INSTALL_DIR = $(PUBLIC)/$(SYSTEM)/bin  # local settings on my machine

# Include path (needs probably no adaption)
INCL    = -I.

# Path for the HTML conversion tools for the default configuration files.
# This path will be compiled into the cc_extract_html program and used for 
# the installation.
# The same variable has to be configured in the cc_manual_to_html script.
# The path must terminate with a slash.
#HTML_DEFAULT_PATH = 
HTML_DEFAULT_PATH = /pub/local/CGAL/Tools/html_config_files/

# Three supported C++ compilers: 
#   Gnu g++ 2.6.3 (or above)
#   Sun Pro CC 4.0.1 (or above)
#   SGI CC 4.0 (or above)
# The compiler call:
#CPP = g++
CPP = CC

# The C++ compiler options (optimization, debugging, warnings, includes)
#CPPOPT = -O -Wall $(INCL)   # Gnu
#CPPOPT = -O +w    $(INCL)   # Sgi
CPPOPT = -O +w    $(INCL)   # Sun

# A C compiler: 
#   Gnu gcc 2.6.3 (or above)
#   System cc
#CC = gcc
CC = cc

# The C compiler options
# CCOPT = $(CPPOPT)   # Gnu
CCOPT = -O $(INCL)  # Sun, SGI

# flex, bison, compile and link commands
FLEX    = flex -t -8
BISON   = bison -d -t -v
COMPILE = $(CPP) $(CPPOPT) -c 
LINK    = $(CPP) $(CPPOPT)

# The flex and bison output might produce a lot of warnings with 
# the C++ compiler (Sgi). Special options switch these messages off.
#XCOMPILE = $(COMPILE)   # Gnu
#XCOMPILE = $(COMPILE) -woff 3203    # SGI
XCOMPILE = $(COMPILE)   # Sun

# Default targets to make:
ALL = cc_extract  cc_build_checker  cc_extract_html  cc_index_sort

# =======================================================
# Installation: Nothing else below this line.
# =======================================================

# all the files that have to be copied into a directory in the $PATH
PROGRAM_FILES = cc_extract  cc_build_checker  cc_extract_html  cc_index_sort

# all the files that have to be moved into a directory in the $PATH
SCRIPT_FILES  = cc_check  cc_manual_to_html

# all the files for the HTML config directory
HTML_CONFIG_FILES = html_config_files/*
                   # cc_anchor_footer cc_anchor_header cc_biblio_footer \
                   # cc_biblio_header cc_index_footer cc_index_header   \
                   # cc_manual_footer cc_manual_header cc_toc_footer    \
                   # cc_toc_header

# default rule
all: $(ALL)


# General Modules
# ==========================
database.o : database.C database.h config.h
	$(COMPILE) -o database.o database.C

# Extraction and Checking
# ==========================
# parser and scanner
# --------------------
extract_lex.c : extract_lex.yy
	$(FLEX) extract_lex.yy >extract_lex.c

extract_syntax.tab.c extract_syntax.tab.h : extract_syntax.y
	$(BISON) extract_syntax.y

# modules
# --------------------
extract_lex.o : extract_lex.c extract_syntax.tab.h database.h config.h
	$(XCOMPILE) -o extract_lex.o extract_lex.c

extract_syntax.o : extract_syntax.tab.c extract_syntax.tab.h \
                   database.h config.h
	$(XCOMPILE) -o extract_syntax.o extract_syntax.tab.c

cc_extract.o : cc_extract.C database.h config.h
	$(COMPILE) -o cc_extract.o cc_extract.C

cc_build_checker.o : cc_build_checker.C database.h config.h
	$(COMPILE) -o cc_build_checker.o cc_build_checker.C

# program
# --------------------
OBJ        = extract_lex.o extract_syntax.o database.o
HEADER     = extract_syntax.tab.h database.h config.h

cc_extract : $(OBJ) $(HEADER) cc_extract.o
	$(LINK) -o cc_extract cc_extract.o $(OBJ)

cc_build_checker : $(OBJ) $(HEADER) cc_build_checker.o
	$(LINK) -o cc_build_checker cc_build_checker.o $(OBJ)


# Converting to HTML
# ==========================
# parser and scanner
# --------------------
html_lex.c : html_lex.yy
	$(FLEX) html_lex.yy >html_lex.c

html_syntax.tab.c html_syntax.tab.h : html_syntax.y
	$(BISON) html_syntax.y

# modules
# --------------------
html_lex.o : html_lex.c html_syntax.tab.h database.h html_config.h
	$(XCOMPILE) -o html_lex.o html_lex.c

html_syntax.o : html_syntax.tab.c html_syntax.tab.h database.h html_config.h
	$(XCOMPILE) -o html_syntax.o html_syntax.tab.c

cc_extract_html.o : cc_extract_html.C database.h html_config.h
	$(COMPILE) -DHTML_DEFAULT_PATH=\"$(HTML_DEFAULT_PATH)\" -o cc_extract_html.o cc_extract_html.C


# program
# --------------------
OBJHTML    = html_lex.o html_syntax.o database.o
HEADERHTML = html_syntax.tab.h database.h html_config.h

cc_extract_html : $(OBJHTML) $(HEADERHTML) cc_extract_html.o
	$(LINK) -o cc_extract_html cc_extract_html.o $(OBJHTML)


# Sort Program to Sort the Index
# ==============================
cc_index_sort : cc_index_sort.c
	$(CC) $(CCOPT) -o cc_index_sort cc_index_sort.c

# Installation Targets
# ====================
install:
	mv $(PROGRAM_FILES) $(INSTALL_DIR)
	cp $(SCRIPT_FILES)  $(INSTALL_DIR)
	cp $(HTML_CONFIG_FILES) $(HTML_DEFAULT_PATH)

# Auxiliary Targets
# =================
clean :
	rm *.o extract_lex.c extract_syntax.output extract_syntax.tab.h \
	       extract_syntax.tab.c html_lex.c html_syntax.output \
	       html_syntax.tab.h html_syntax.tab.c

cleanall : clean
	rm cc_build_checker cc_extract  cc_extract_html
	if (-d Templates.DB) rm -r Templates.DB
	if (-d ii_files) rm -r ii_files
