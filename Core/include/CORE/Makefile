# file: Makefile (for inc)
#	CORE Library version 1.4
#	$Id$

# =============================================================
# PATHS =======================================================
# =============================================================

CORE_PATH=..

INCLUDE=-I${CORE_PATH}/inc

LIBPATH= 

save:
	-mkdir RCS; ci -l -m"make save" *.h Makefile*

clean:

veryclean: clean

REPLACE_STRING="s/1.6, June 2003/1.7, August 2004/;s/1995-2003/1995-2004/"

chver:
	@sed -e ${REPLACE_STRING} ../CORE.h > x.tmp; mv x.tmp ../CORE.h
	@for f in *.h; do \
	  sed -e ${REPLACE_STRING} $$f > x.tmp; mv x.tmp $$f; \
	done
	@for f in poly/*.h; do \
	  sed -e ${REPLACE_STRING} $$f > x.tmp; mv x.tmp $$f; \
	done
	@for f in poly/*.tcc; do \
	  sed -e ${REPLACE_STRING} $$f > x.tmp; mv x.tmp $$f; \
	done
	@for f in geom2d/*.h; do \
	  sed -e ${REPLACE_STRING} $$f > x.tmp; mv x.tmp $$f; \
	done
	@for f in geom3d/*.h; do \
	  sed -e ${REPLACE_STRING} $$f > x.tmp; mv x.tmp $$f; \
	done

