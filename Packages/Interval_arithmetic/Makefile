# Makefile for a CGAL package.

.PHONY: doc specific test submit Version

Package:=$(shell basename `cat CVS/Repository`)
TAG_VERSION:=v$(shell echo $(VERSION) | sed -e "s/\./_/g")
CGAL_VERSION:="$(VERSION) ($(shell date '+%e %b %Y'))"
WWW_base=/net/servers/www-sop/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

all:
	@echo "-----------------------------------"
	@echo "---- edit changes.txt          ----"
	@echo "---- make VERSION=X.XX submit  ----"
	@echo "-----------------------------------"

Version:
	echo $(CGAL_VERSION) > version
	cvs commit -m "- Bump for $(CGAL_VERSION)" version
	cvs tag -c $(TAG_VERSION)

submit: specific clean Version $(Package).tar.gz
	cp $(Package).tar.gz doc_ps/$(Package).ps version description.txt \
	   changes.txt $(WWW_dir)
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl $(USER) < $(TMP_file)
	rm -f $(TMP_file)

test:
	$(MAKE) -C test/$(Package)

clean:
	$(MAKE) -C doc_tex/support/NumberTypeSupport clean
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz src/.pure src/*.o

$(Package).tar.gz: doc
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile --exclude=TODO \
	--exclude="*".swp --exclude=tst_Static_filter_error.C \
	--exclude=tst_static_adaptatif.C \
	version description.txt changes.txt maintainer \
	src include doc_tex doc_ps test scripts

doc: doc_ps/$(Package).ps

doc_ps/$(Package).ps: doc_tex/support/NumberTypeSupport/Interval.tex
	$(MAKE) -C doc_tex/support/NumberTypeSupport
	cp doc_tex/support/NumberTypeSupport/$(Package).ps doc_ps/

# The following is specific to the IA package.

specific:
	$(MAKE) -C include/CGAL/Arithmetic_filter
