#
# This is the Makefile responsible for constructing the filtered predicates
# files from the generic ones in $CGAL/include/CGAL.
#
# Maybe you need GNU Make, I'm not sure.
#
# Sylvain Pion, 1998.

# Tell where the script is.
SCRIPT = ../../../filtered_predicate_converter

# Tell which files must be parsed.
OBJECTS = template_predicates_on_ft.h

#### Don't touch below.

% :: ../% $(SCRIPT)
	$(SCRIPT) < $< > $@

all: $(OBJECTS)

clean:
	rm -f $(OBJECTS)
