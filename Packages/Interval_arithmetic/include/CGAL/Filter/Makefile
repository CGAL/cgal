#
# This is the Makefile responsible for constructing the filtered predicates
# files from the generic ones in $CGAL/include/CGAL.
#
# Maybe you need GNU Make, I don't know.
#
# Sylvain Pion, 1998.

# Tell where the script is.
SCRIPT = ../../../script/filtered_predicate_converter

# Tell which files must be parsed.
OBJECTS = predicates_on_ftC2.h

#### Don't touch below.

$(OBJECTS): %.h: $(CGAL)/include/CGAL/%.h $(SCRIPT)
	$(SCRIPT) < $< > $@

all: $(OBJECTS)

clean:
	rm -f $(OBJECTS)
