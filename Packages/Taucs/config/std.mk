DIRSRC=src$(PATHSEP)

DIRPROGS=progs$(PATHSEP)

DIREXTSRC=external$(PATHSEP)src$(PATHSEP)
DIREXTLIB=external$(PATHSEP)lib$(PATHSEP)

DIRLIB=lib$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)
DIROBJ=obj$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)
DIREXE=bin$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)
DIRBLD=build$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)

STDINCS = -I $(DIRSRC) -I $(DIRBLD) -I $(DIREXTSRC) $(OSTYPE_INCS) $(OSTYPE_VARIANT_INCS)

STDDEFS = \
  $(OSTYPE_DEFS) \
  $(OSTYPE_VARIANT_DEFS) \
  -DMACHTYPE_$(MACHTYPE)

STDDEPS = $(DIROBJ)exists.log \
          $(DIREXE)exists.log \
          $(DIRLIB)exists.log

# These are for nmake, which refuses to make
# a target (exists.log) if the directory where
# is it to reside does not exist. So these
# tell nmake how to make the directories.
obj:
	- mkdir obj
bin:
	- mkdir bin
lib:
	- mkdir lib
$(DIROBJ): obj
	- mkdir $(DIROBJ)
$(DIREXE): bin
	- mkdir $(DIREXE)
$(DIRLIB): lib
	- mkdir $(DIRLIB)

$(DIROBJ)exists.log:
	- mkdir obj
	- mkdir $(DIROBJ)
	- echo exists > $(DIROBJ)exists.log

$(DIREXE)exists.log:
	- mkdir bin
	- mkdir $(DIREXE)
	- echo exists > $(DIREXE)exists.log

$(DIRLIB)exists.log:
	- mkdir lib
	- mkdir $(DIRLIB)
	- echo exists > $(DIRLIB)exists.log

clean:
	- $(RM) lib$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)*
	- $(RM) obj$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)*
	- $(RM) bin$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)*
	- $(RM) lib$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)*.*
	- $(RM) obj$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)*.*
	- $(RM) bin$(PATHSEP)$(OSTYPE)$(OSTYPE_VARIANT)$(PATHSEP)*.*
	- $(RM) $(DIRBLD)taucs_config_tests.h
	- rmdir $(DIROBJ)
	- rmdir $(DIREXE)
	- rmdir $(DIRLIB)

distclean:
	- $(RM) *.tar *.zip
	- $(RM) configurator$(PATHSEP)configurator
	- $(RM) configurator$(PATHSEP)configurator.exe
	- $(RM) doc$(PATHSEP)taucs.tex
	- $(RM) doc$(PATHSEP)taucs.dvi
	- $(RM) doc$(PATHSEP)taucs.log
	- $(RM) doc$(PATHSEP)taucs.aux
	- $(RM) doc$(PATHSEP)taucs.bbl
	- $(RM) doc$(PATHSEP)taucs.blg
	- $(RM) doc$(PATHSEP)taucs.ps
	- $(RM) doc$(PATHSEP)taucs.toc
	- $(RM) .lastconf
	- $(RM) testscript.log
	- $(RM) lib
	- $(RM) bin
	- $(RM) obj
	- $(RM) build

FILES = \
	config$(PATHSEP)*.mk \
	configurator$(PATHSEP)makefile.* \
	configurator$(PATHSEP)taucs_*.* \
	doc$(PATHSEP)taucs.pdf \
	doc$(PATHSEP)cilk-mf.pdf \
	$(DIREXTSRC)*.h \
	$(DIREXTSRC)*.c \
	$(DIREXTSRC)*.f \
	matlab$(PATHSEP)*.m \
	$(DIRPROGS)*.c \
	$(DIRPROGS)mktestscript.sh \
	testscript \
	testscript.bat \
	$(DIRSRC)*.c \
	$(DIRSRC)*.h \
        configure.bat configure 

LIBFILES = $(DIREXTLIB)*

testscript: $(DIRPROGS)mktestscript.sh
	/bin/sh $(DIRPROGS)mktestscript.sh

testscript.bat: $(DIRPROGS)mktestscript.sh
	/bin/sh $(DIRPROGS)mktestscript.sh

tar: $(FILES) $(LIBFILES)
	@echo "Creating distribution files"
	-$(RM) taucs.tgz taucs.zip taucs_full.tgz taucs_full.zip
	tar zcpf taucs.tgz      $(FILES)
	tar zcpf taucs_full.tgz $(FILES) $(LIBFILES)
	zip -r   taucs.zip      $(FILES)
	zip -r   taucs_full.zip $(FILES) $(LIBFILES)
	chmod 644 *.zip *.tgz

#	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
#	-DTAUCS_BLAS_UNDERSCORE \
#	$(COUTFLG)$(DIROBJ)taucs_blas_link_test_yes$(OBJEXT) $(DIRPROGS)taucs_blas_link_test.c
#	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
#	$(COUTFLG)$(DIROBJ)taucs_blas_link_test_no$(OBJEXT) $(DIRPROGS)taucs_blas_link_test.c
#	- $(LD) $(LDFLAGS) \
#	$(LOUTFLG)$(DIREXE)taucs_blas_link_test_yes$(EXEEXT) \
#	$(DIROBJ)taucs_blas_link_test_yes$(OBJEXT) \
#	$(LIBS)
#	- $(LD) $(LDFLAGS) \
#	$(LOUTFLG)$(DIREXE)taucs_blas_link_test_no$(EXEEXT) \
#	$(DIROBJ)taucs_blas_link_test_no$(OBJEXT) \
#	$(LIBS)
#	- $(DIREXE)taucs_blas_link_test_yes$(EXEEXT) $(DIRBLD)taucs_config_tests.h
#	- $(DIREXE)taucs_blas_link_test_no$(EXEEXT) $(DIRBLD)taucs_config_tests.h

$(DIRBLD)OLD_taucs_config_tests.h: $(DIRPROGS)taucs_blas_link_test.c $(DIREXE)exists.log $(DIROBJ)exists.log
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	$(COUTFLG)$(DIROBJ)taucs_blas_underscore_test$(OBJEXT) $(DIRPROGS)taucs_blas_underscore_test.c
	- $(LD) $(LDFLAGS) \
	$(LOUTFLG)$(DIREXE)taucs_blas_underscore_test$(EXEEXT) \
	$(DIROBJ)taucs_blas_underscore_test$(OBJEXT) \
	$(LIBS)
	- $(DIREXE)taucs_blas_underscore_test$(EXEEXT) $(DIRBLD)taucs_config_tests.h
	$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	$(COUTFLG)$(DIROBJ)taucs_blas_nounderscore_test$(OBJEXT) $(DIRPROGS)taucs_blas_nounderscore_test.c
	- $(LD) $(LDFLAGS) \
	$(LOUTFLG)$(DIREXE)taucs_blas_nounderscore_test$(EXEEXT) \
	$(DIROBJ)taucs_blas_nounderscore_test$(OBJEXT) \
	$(LIBS)
	- $(DIREXE)taucs_blas_nounderscore_test$(EXEEXT) $(DIRBLD)taucs_config_tests.h
	- $(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \
	$(COUTFLG)$(DIROBJ)taucs_c99_complex_test$(OBJEXT) $(DIRPROGS)taucs_c99_complex_test.c
	- $(LD) $(LDFLAGS) \
	$(LOUTFLG)$(DIREXE)taucs_c99_complex_test$(EXEEXT) \
	$(DIROBJ)taucs_c99_complex_test$(OBJEXT) $(LIBS)
	- $(DIREXE)taucs_c99_complex_test$(EXEEXT) $(DIRBLD)taucs_config_tests.h
	- $(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \
	$(CILKOUTFLG)$(DIROBJ)taucs_cilk_test$(OBJEXT) $(DIRPROGS)taucs_cilk_test.c
	- $(LD) $(LDFLAGS) \
	$(LOUTFLG)$(DIREXE)taucs_cilk_test$(EXEEXT) \
	$(DIROBJ)taucs_cilk_test$(OBJEXT) $(LIBS)
	- $(DIREXE)taucs_cilk_test$(EXEEXT) $(DIRBLD)taucs_config_tests.h

STD_PRE_TARGETS=



