CGAL_URL := $(shell ./get_cgal_trunk_url)

.PHONY: \
        SOURCES/Manual_tools.tar.gz \
        SOURCES/Manual.tar.gz

prepare-CGAL-manual-tools-sources: SOURCES/Manual_tools.tar.gz \
                                   SOURCES/Manual.tar.gz

SOURCES/Manual_tools.tar.gz:
	cd SOURCES && \
	  svn export --force $(CGAL_URL)/Manual_tools && \
	  tar czf Manual_tools.tar.gz Manual_tools/
	rm -rf SOURCES/Manual_tools/

SOURCES/Manual.tar.gz:
	cd SOURCES && \
	  svn export --force $(CGAL_URL)/Manual && \
	  tar czf Manual.tar.gz -C Manual/ .
	rm -rf SOURCES/Manual/

-include my_makefile
