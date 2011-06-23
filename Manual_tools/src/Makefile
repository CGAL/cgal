#  Makefile
#  ===================================
#     V 1.0    17.08.1995   Lutz Kettner
#              Project: CGAL merger tool for the specification task
#     V 1.1    22.11.1996   Lutz Kettner
#              Project: tools around the CGAL manual writing task.
#              extract, check, HTML conversion.
#     V 1.2    03.12.1996   Lutz Kettner
#              Tools renamed from CGAL to CC. Relates to cc_manual.sty R2.1.
#     V 1.3    03.08.1998   Lutz Kettner
#              Simplified installation. LaTeX converter.

# =====================================================================
# Begin of Installation Section: ...
# =====================================================================
# Choose appropriate values in the environment or uncomment and
# set the values here. See the INSTALLATION file for more details.

# Path to the configuration files of the latex_converter:
# Used during installation to copy config files and compiled
# into the cc_extract_html program as default config path.
# The same variable is used in the latex_to_html script.
# LATEX_CONV_CONFIG = ${PUBLIC}/CGAL/Tools/latex_converter_config

# Path where to copy/move the executables during installation:
# LATEX_CONV_BIN = ${PUBLIC}/${SYSTEM}/bin

# Analog of the TEXINPUTS environment variable for the cc_extract_html
# program. It should contain a colon-separated list of directories where
# the input files can be found. If this is not set, the default is that
# files are looked for in the current directory. If this is set, then
# directories are searched in the order provided in this variable (so "."
# should be included in the list if some input files are in the current
# directory).
LATEX_CONV_INPUTS = .

# C++ Compiler: tested compilers are: Gnu g++ 3.3 and 3.4
# CXX = g++-3.4

# C++ Compiler Options:
# CXXFLAGS = -v
CXXFLAGS = -g -O0

# C Compiler: Gnu gcc 2.6.3 (or above) or the usual system cc:
# CC = gcc

# C Preprocessor Flags: (e.g. -g or -O or STL include path)
# CPPFLAGS = 
CPPFLAGS = 


# Linker Flags:
# LDFLAGS = 

# Flex and Bison Tools (not supposed to be set in the environment):
FLEX    = flex -t -8 -Cfe
BISON   = bison -d -t -v


# =====================================================================
# ... End of Installation Section.
# =====================================================================

# =====================================================================
# File Lists and Default Target:
# =====================================================================

# Files that will be moved into LATEX_CONV_BIN upon installation:
PROGRAM_FILES =		cc_extract_html			\
			cc_index_link                   \
	        	cc_extract_images		\
			cc_anchor_filter                \
                        cc_remove_unwanted_links


# Files that will be copied into LATEX_CONV_BIN upon installation:
SCRIPT_FILES  = 	latex_to_html			\
			cc_extract_include		

# Files that will be copied into LATEX_CONV_CONFIG upon installation:
HTML_STYLE_IMAGES =	latex_converter_config/html/style/*.gif
HTML_CONFIG_FILES =	latex_converter_config/html/*.sty	\
		 	latex_converter_config/html/*.css	\
			latex_converter_config/html/*.mst	
BIB_CONFIG_FILES =	latex_converter_config/html/Biblio/how_to_cite_cgal.html \
			latex_converter_config/html/Biblio/how_to_cite_cgal.bib \
			latex_converter_config/html/Biblio/cgal.bib
GIF_CONFIG_FILES  =	latex_converter_config/gif/*.gif
REF_CONFIG_FILES  =	latex_converter_config/ref_pages/*.tex


# Default targets to make:
ALL =   ${PROGRAM_FILES}

# default rule
all: $(ALL) 


# =====================================================================
# Dependancies and Specific Production Rules:
# =====================================================================

# Extraction and Checking
# ==========================

# modules
# --------------------

basic.o : basic.h

buffer.o : buffer.cpp buffer.h config.h basic.h


# LaTeX Converter
# ==========================

# modules
# --------------------
include Makefile_depend

cc_extract_html.o : cc_extract_html.cpp
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} -c cc_extract_html.cpp

cc_anchor_filter.o : cc_anchor_filter.cpp
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} -c cc_anchor_filter.cpp

cc_remove_unwanted_links.o : cc_remove_unwanted_links.cpp
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} -c cc_remove_unwanted_links.cpp

# program
# --------------------
OBJHTML    = lexer.o    \
             syntax.o \
             basic.o             \
             buffer.o            \
             macro_dictionary.o  \
             internal_macros.o   \
             error.o        \
             string_conversion.o \
             cpp_formatting.o    \
             input.o             \
             output.o            \
             cc_extract_html.o	\
	     regex++.o		\
	     regex_dictionary.o

OBJIMAGES = cc_extract_images_lexer.cpp cc_extract_images.o


# configure input.cpp with the right parameters
input.o: input.cpp
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} \
		-DLATEX_CONVERTER_CONFIG=\"${LATEX_CONV_CONFIG}\" \
		-DLATEX_CONVERTER_INPUTS=\"${LATEX_CONV_INPUTS}\" \
		-c input.cpp


cc_extract_html : $(OBJHTML)
	${CXX} ${CXXFLAGS} -o $@ ${OBJHTML} ${LDFLAGS}

cc_extract_images : $(OBJIMAGES)
	${CXX} ${CXXFLAGS} -o $@ ${OBJIMAGES} ${LDFLAGS}

cc_index_link : cc_index_link.o
	${CXX} ${CXXFLAGS} -o $@ $@.o ${LD_FLAGS}

cc_anchor_filter : cc_anchor_filter.o
	${CXX} ${CXXFLAGS} -o $@ $@.o ${LD_FLAGS}

cc_remove_unwanted_links : cc_remove_unwanted_links.o
	${CXX} ${CXXFLAGS} -o $@ $@.o ${LD_FLAGS}

# =====================================================================
# Generic Production Rules:
# =====================================================================

%.cpp: %.lex
	${FLEX} $*.lex > $*.cpp

%.cpp %.tab.h: %.y
	${BISON} $*.y
	mv $*.tab.c $*.cpp

%.o: %.cpp
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} -c $*.cpp

%.o: %.c
	${CC} -I. ${CPPFLAGS} -c $*.c

%: %.c
	${CC} -I. ${CPPFLAGS} -o $* $*.c ${LDFLAGS}

# =====================================================================
# Installation Targets:
# =====================================================================
rminstall:
	test -d "${LATEX_CONV_CONFIG}" || rm -fr "${LATEX_CONV_CONFIG}"

pathinstall:
	test -d "${LATEX_CONV_CONFIG}"      || mkdir "${LATEX_CONV_CONFIG}"
	test -d "${LATEX_CONV_CONFIG}/html" || mkdir "${LATEX_CONV_CONFIG}/html"
	test -d "${LATEX_CONV_CONFIG}/html/style" || mkdir "${LATEX_CONV_CONFIG}/html/style"
	test -d "${LATEX_CONV_CONFIG}/html/Biblio" || mkdir "${LATEX_CONV_CONFIG}/html/Biblio"
	test -d "${LATEX_CONV_CONFIG}/gif"  || mkdir "${LATEX_CONV_CONFIG}/gif"
	test -d "${LATEX_CONV_CONFIG}/ref_pages" || mkdir "${LATEX_CONV_CONFIG}/ref_pages"

gifinstall:
	install -m 644 ${GIF_CONFIG_FILES}  "${LATEX_CONV_CONFIG}/gif"

refinstall:
	install -m 644 ${REF_CONFIG_FILES}  "${LATEX_CONV_CONFIG}/ref_pages"

cpinstall:
	install ${SCRIPT_FILES}  "${LATEX_CONV_BIN}"
	install -m 644 ${HTML_CONFIG_FILES} "${LATEX_CONV_CONFIG}/html"
	install -m 644 ${HTML_STYLE_IMAGES} "${LATEX_CONV_CONFIG}/html/style"
	install -m 644 ${BIB_CONFIG_FILES} "${LATEX_CONV_CONFIG}/html/Biblio"

mvinstall: ${PROGRAM_FILES}
	install ${PROGRAM_FILES} ${LATEX_CONV_BIN}

install: rminstall pathinstall cpinstall mvinstall refinstall gifinstall


# =====================================================================
# Target to Create Dependencies Automatically:
# =====================================================================

DEPENDSRC  = lexer.cpp    \
             syntax.cpp \
             basic.cpp             \
             buffer.cpp            \
             macro_dictionary.cpp  \
             internal_macros.cpp   \
             error.cpp        \
             string_conversion.cpp \
             cpp_formatting.cpp    \
             input.cpp             \
             output.cpp            \
             cc_extract_html.cpp


# Prepared for SunPro CC. Change -xM1 to -M, -MM, or -MD, for Gnu g++.

#depend: ${DEPENDSRC}
#	# Works only with SunPro CC
#	CC -I. ${CPPFLAGS} ${CXXFLAGS} -xM1 ${DEPENDSRC} | grep -v "STL/SUN" | grep -v "opt/SUNWspro" | sort | uniq > Makefile_depend

depend: ${DEPENDSRC}
	${CXX} -MM -I. ${CPPFLAGS} ${CXXFLAGS} ${DEPENDSRC} > Makefile_depend

# =====================================================================
# Auxiliary Targets:
# =====================================================================
clean :
	-rm -f	*.o *~ *.bak	\
		lexer.cpp		\
		syntax.output	\
		syntax.tab.h	\
		syntax.cpp	\
		cc_remove_unwanted_links.cpp \
		cc_anchor_filter.cpp \
		cc_extract_images_lexer.cpp

cleanall : clean
	-rm -f  $(PROGRAM_FILES)
	-rm -fr Templates.DB
	-rm -fr ii_files
