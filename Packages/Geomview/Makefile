# Makefile for the Geomwiew package.

.PHONY: demo doc submit test Version

Package:=$(shell basename `cat CVS/Repository`)
TAG_VERSION:=v$(shell echo $(VERSION) | sed -e "s/\./_/g")
CGAL_VERSION:="$(VERSION) ($(shell date '+%e %b %Y'))"
WWW_CGAL=/net/servers/www-sop/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/$(Package)
TMP_file=/tmp/submit-mail

all:
	@echo "----------------------------------"
	@echo "---- edit changes.txt ------------"
	@echo "---- make VERSION=X.Y submit -----"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version
	cvs commit -m "- Bump for $(CGAL_VERSION)" version
	cvs tag -c $(TAG_VERSION)

$(Package).tar.gz: doc
	$(MAKE) clean
	tar -zcvf $(Package).tar.gz --exclude="*".swp \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	version maintainer description.txt changes.txt \
	demo include src test doc_tex doc_ps

submit: Version $(Package).tar.gz
	cp -f $(Package).tar.gz version doc_ps/$(Package).ps.gz *.txt $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Geomview/Geomview.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl $(USER) < $(TMP_file)
	rm -f $(TMP_file)

doc:
	$(MAKE) -C doc_tex/support/Geomview all
	cp -f doc_tex/support/Geomview/Geomview.ps.gz doc_ps

test:
	$(MAKE) -C test/Geomview all

demo:
	$(MAKE) -C demo/Geomview all

clean:
	$(MAKE) -C test/Geomview clean
	$(MAKE) -C demo/Geomview clean
	$(MAKE) -C doc_tex/support/Geomview clean
