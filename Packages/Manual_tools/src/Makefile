#  Makefile
#  ===================================

nodefault:
	# No default targets.
	# Use the system specific compiler makefiles:
	#     Makefile.sun, Makefile.sgi, Makefile.gcc
	# Possible other targets:
	#     install, cpinstall, clean, cleanall


# Installation directory (where to put the executables)
# Only necessary if the make target 'install' is used for installation.
#INSTALL_DIR = /usr/bin
INSTALL_DIR = $(PUBLIC)/$(SYSTEM)/bin  # local settings on my machine

# =======================================================

# all the files that have to be copied into a directory in the $PATH
PROGRAM_FILES = cc_extract  cc_build_checker  cc_extract_html  cc_index_sort

# all the files that have to be moved into a directory in the $PATH
SCRIPT_FILES  = cc_check  cc_manual_to_html  cc_extract_include \
                cc_patch_anchor_pages  cc_patch_anchor_toc_index \
                cc_patch_anchor_filter

# all the files for the HTML config directory
HTML_CONFIG_FILES = html_config_files/cc_*


# Installation Targets
# ====================
install:
	mv -f $(PROGRAM_FILES) $(INSTALL_DIR)
	cp $(SCRIPT_FILES)  $(INSTALL_DIR)
	-cp $(HTML_CONFIG_FILES) $(HTML_DEFAULT_PATH) 2>/dev/null

cpinstall:
	cp $(SCRIPT_FILES)  $(INSTALL_DIR)
	-cp $(HTML_CONFIG_FILES) $(HTML_DEFAULT_PATH) 2>/dev/null

# Auxiliary Targets
# =================
clean :
	rm -f *.o extract_lex.c extract_syntax.output extract_syntax.tab.h \
	       extract_syntax.tab.c html_lex.c html_syntax.output \
	       html_syntax.tab.h html_syntax.tab.c

cleanall : clean
	rm -f cc_build_checker cc_extract  cc_extract_html
	if (-d Templates.DB) rm -r Templates.DB
	if (-d ii_files) rm -r ii_files
