#  Makefile for CGAL's CLN package. [Sylvain Pion]

.PHONY: doc package test submit version tag

Package:=$(shell basename `cat CVS/Repository`)
TAG_VERSION:=v$(shell echo $(VERSION) | sed -e "s/\./_/g")
CGAL_VERSION:="$(VERSION) ($(shell date '+%e %b %Y'))"
WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: doc clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make test                ----"
	@echo "---- make VERSION=1.0 version ----"
	@echo "---- cvs commit               ----"
	@echo "---- make VERSION=1.0 tag     ----"
	@echo "---- make package             ----"
	@echo "---- make submit              ----"
	@echo "----------------------------------"

version:
	echo $(CGAL_VERSION) > version

tag:
	cvs tag -c $(TAG_VERSION)

submit: package
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

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile --exclude=TODO \
	--exclude="*".swp \
	version description.txt changes.txt maintainer \
	doc_tex doc_ps include examples src test

doc:
	$(MAKE) -C doc_tex/support/NumberTypeSupport
