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
# The same variable is used in the cc_manual_to_html script.
# LATEX_CONV_CONFIG = ${PUBLIC}/CGAL/Tools/latex_converter_config

# Path where to copy/move the executables during installation:
# LATEX_CONV_BIN = ${PUBLIC}/${SYSTEM}/bin


# C++ Compiler: tested compilers are: Gnu g++ 2.8.1, SGI CC 7.2.
# CXX = g++

# C++ Compiler Options:
# CXXFLAGS = -v

# C Compiler: Gnu gcc 2.6.3 (or above) or the usual system cc:
# CC = gcc

# C Preprocessor Flags: (e.g. -g or -O or STL include path)
# CPPFLAGS = -g

# Linker Flags:
# LDFLAGS = 

# Flex and Bison Tools (not supposed to be set in the environment):
FLEX    = flex -t -8
BISON   = bison -d -t -v


# =====================================================================
# ... End of Installation Section.
# =====================================================================

# =====================================================================
# File Lists and Default Target:
# =====================================================================

# Files that will be moved into LATEX_CONV_BIN upon installation:
PROGRAM_FILES =		cc_extract			\
			cc_build_checker		\
			cc_extract_html			\
	        	cc_extract_images		\
			cc_index_sort
#			cc_unique

# Files that will be copied into LATEX_CONV_BIN upon installation:
SCRIPT_FILES  = 	cc_check			\
			cc_manual_to_html		\
			cc_extract_include		\
			cc_patch_anchor_pages		\
			cc_patch_anchor_toc_index	\
			cc_patch_anchor_filter          \
                        cc_remove_unwanted_links

# Files that will be copied into LATEX_CONV_CONFIG upon installation:
HTML_CONFIG_FILES =	latex_converter_config/html/*.sty	\
			latex_converter_config/html/cc_*
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

buffer.o : buffer.C buffer.h config.h basic.h

extract_lex.o : extract_syntax.tab.h buffer.h config.h basic.h

extract_syntax.o : extract_syntax.tab.h buffer.h config.h basic.h

cc_extract.o : buffer.h config.h basic.h

cc_build_checker.o : buffer.h config.h basic.h

# program
# --------------------
OBJ        = extract_lex.o extract_syntax.o buffer.o basic.o
HEADER     = extract_syntax.tab.h buffer.h config.h basic.h

cc_extract : ${OBJ} ${HEADER} cc_extract.o
	${CXX} ${CXXFLAGS} -o $@ $@.o ${OBJ} ${LDFLAGS}

cc_build_checker : ${OBJ} ${HEADER} cc_build_checker.o
	${CXX} ${CXXFLAGS} -o $@ $@.o ${OBJ} ${LDFLAGS}

texinputs : texinputs.o
	${CXX} ${CXXFLAGS} -o $@ $@.o ${LDFLAGS}

# LaTeX Converter
# ==========================

# modules
# --------------------
include Makefile_depend

cc_extract_html.o : cc_extract_html.C
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} \
	-DLATEX_CONVERTER_CONFIG=\"${LATEX_CONV_CONFIG}\" \
	-c cc_extract_html.C

# program
# --------------------
OBJHTML    = html_lex.o          \
             html_syntax.o       \
             basic.o             \
             buffer.o            \
             macro_dictionary.o  \
             internal_macros.o   \
             html_error.o        \
             string_conversion.o \
             cpp_formatting.o    \
             output.o            \
             cc_extract_html.o

cc_extract_html : $(OBJHTML)
	${CXX} ${CXXFLAGS} -o $@ ${OBJHTML} ${LDFLAGS}


# =====================================================================
# Generic Production Rules:
# =====================================================================

%.C: %.yy
	${FLEX} $*.yy > $*.C

%.C %.tab.h: %.y
	${BISON} $*.y
	mv $*.tab.c $*.C

%.o: %.C
	${CXX} -I. ${CPPFLAGS} ${CXXFLAGS} -c $*.C

%.o: %.c
	${CC} -I. ${CPPFLAGS} -c $*.c

%: %.c
	${CC} -I. ${CPPFLAGS} -o $* $*.c ${LDFLAGS}

# =====================================================================
# Installation Targets:
# =====================================================================
rminstall:
	-rm -f ${LATEX_CONV_CONFIG}/*
	-rm -f ${LATEX_CONV_CONFIG}/html/*
	-rm -f ${LATEX_CONV_CONFIG}/gif/*
	-rm -f ${LATEX_CONV_CONFIG}/ref_pages/*

pathinstall:
	-mkdir ${LATEX_CONV_CONFIG}
	-mkdir ${LATEX_CONV_CONFIG}/html
	-mkdir ${LATEX_CONV_CONFIG}/gif
	-mkdir ${LATEX_CONV_CONFIG}/ref_pages

gifinstall:
	-cp ${GIF_CONFIG_FILES}  ${LATEX_CONV_CONFIG}/gif

refinstall:
	-cp ${REF_CONFIG_FILES}  ${LATEX_CONV_CONFIG}/ref_pages

cpinstall:
	cp ${SCRIPT_FILES}  ${LATEX_CONV_BIN}
	-cp ${HTML_CONFIG_FILES} ${LATEX_CONV_CONFIG}/html

mvinstall:
	mv -f ${PROGRAM_FILES} ${LATEX_CONV_BIN}

install: rminstall pathinstall cpinstall mvinstall refinstall gifinstall


# =====================================================================
# Target to Create Dependencies Automatically:
# =====================================================================

DEPENDSRC  = html_lex.C          \
             html_syntax.C       \
             basic.C             \
             buffer.C            \
             macro_dictionary.C  \
             internal_macros.C   \
             html_error.C        \
             string_conversion.C \
             cpp_formatting.C    \
             output.C            \
             cc_extract_html.C


# Prepared for SunPro CC. Change -xM1 to -M, -MM, or -MD, for Gnu g++.

depend: ${DEPENDSRC}
	# Works only with SunPro CC
	CC -I. ${CPPFLAGS} ${CXXFLAGS} -xM1 ${DEPENDSRC} | grep -v "STL/SUN" | grep -v "opt/SUNWspro" | sort | uniq > Makefile_depend


# =====================================================================
# Auxiliary Targets:
# =====================================================================
clean :
	-rm -f	*.o			\
		extract_lex.C		\
		extract_syntax.output	\
		extract_syntax.tab.h	\
		extract_syntax.C	\
		html_lex.C		\
		html_syntax.output	\
		html_syntax.tab.h	\
		html_syntax.C

cleanall : clean
	-rm -f  $(PROGRAM_FILES)
	-rm -fr Templates.DB
	-rm -fr ii_files
