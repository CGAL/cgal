# Makefile for a CGAL package.

.PHONY: submit specific

TAG_VERSION:=v$(shell echo $(VERSION) | sed -e "s/\./_/g")
CGAL_VERSION:='$(VERSION) ($(shell date '+%e %b %Y'))'

all:
	@echo "-----------------------------------"
	@echo "---- edit changes.txt          ----"
	@echo -n "---- the current version file contains : "
	@cat version
	@echo "---- make VERSION=X.XX submit  ----"
	@echo "-----------------------------------"

submit: specific
	echo $(CGAL_VERSION) > version
	cvs commit -m "Bump for $(CGAL_VERSION)" version
	cvs tag -c $(TAG_VERSION)
	cvs tag -F -c current_submission

specific:
	$(MAKE) -C include/CGAL/Arithmetic_filter
