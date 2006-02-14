#
# This is the Makefile responsible for constructing the filtered predicates
# files from the generic ones in $CGAL/include/CGAL.
#
# Maybe you need GNU Make, I don't know.
#
# Sylvain Pion, 1998-1999.

# Tell where the script is.
SCRIPT = ../../../scripts/filtered_predicates_generator.pl
BASE = $(CGAL)/include/CGAL

# Tell which files must be treated.
OBJECTS = predicates/kernel_ftC2.h \
          predicates/kernel_ftC3.h \
          predicates/sign_of_determinant.h \
          predicates/Regular_triangulation_ftC2.h \
          predicates/Regular_triangulation_rtH2.h \
          predicates/Regular_triangulation_ftC3.h \
          predicates/Regular_triangulation_rtH3.h \
          predicates/in_smallest_orthogonalcircle_ftC2.h

ALL_DEP = $(BASE)/predicates/sign_of_determinant.h

#### Don't touch below.

% :: $(BASE)/% $(SCRIPT)
	perl $(SCRIPT) -p -d $(ALL_DEP) -l static_infos/$@ -i $< -o $@

predicates/% : $(BASE)/predicates/% $(SCRIPT)
	perl $(SCRIPT) -p -d $(ALL_DEP) -l static_infos/$@ -i $< -o $@

all: $(OBJECTS)

clean:
	rm -f $(OBJECTS)
