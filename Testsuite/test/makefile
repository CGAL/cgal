# A makefile which calls run_testsuite over all directories.
#
# Author : Sylvain Pion

# The "all" target is split in several parts (based on the first
# character of each directory) to allow parallelization using "make -j".

all: all0 all1 all2 all3 all4 all5 all6 all7 all8

all0:
	@+for DIR in [0-9_aAbBcC]* ; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all1:
	@+for DIR in [dDeEfF]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all2:
	@+for DIR in [gGhHiI]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all3:
	@+for DIR in [jJkKlL]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all4:
	@+for DIR in [mMnNoO]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all5:
	@+for DIR in [pPqQrR]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all6:
	@+for DIR in [sStTuU]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all7:
	@+for DIR in [vVwWxX]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

all8:
	@+for DIR in [yYzZ]*; do \
	  if test -d "$$DIR"; then \
	    ( ./run_testsuite "$$DIR" ) \
	  fi; \
	done

clean:
	@for DIR in *; do \
	  if test -d "$$DIR"; then \
	    ( cd "$$DIR" && $(MAKE) clean ) \
	  fi; \
	done
