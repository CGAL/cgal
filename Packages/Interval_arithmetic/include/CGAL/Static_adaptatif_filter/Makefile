#
# This is the Makefile responsible for constructing the filtered predicates
# files from the generic ones in $CGAL/include/CGAL.
#
# Maybe you need GNU Make, I don't know.
#
# Sylvain Pion, 1998-1999.

# Tell where the script is.
SCRIPT = ../../../examples/Interval_arithmetic/filtered_predicate_converter
BASE = $(CGAL)/include/CGAL

# Tell which files must be treated.
OBJECTS = predicates_on_ftC2.h \
          predicates_on_ftC3.h \
          predicates/sign_of_determinant.h \
          predicates/Regular_triangulation_ftC2.h \
          predicates/Regular_triangulation_rtH2.h \
          predicates/in_smallest_orthogonalcircle_ftC2.h

#### Don't touch below.

% :: $(BASE)/% $(SCRIPT)
	$(SCRIPT) -d $(BASE)/predicates/sign_of_determinant.h -p -s < $< > $@

predicates/% : $(BASE)/predicates/% $(SCRIPT)
	$(SCRIPT) -d $(BASE)/predicates/sign_of_determinant.h -p -s < $< > $@

all: $(OBJECTS)

clean:
	rm -f $(OBJECTS)
