#  Makefile for CGAL's Triangulation3 package.

Package=Triangulation3
CVS_MODULE_NAME=$(Package)

.PHONY: doc package test submit version rtag

WWW_base=/net/servers/www-sop/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: doc clean $(Package).tar.gz
	@echo "-----------------------------------"
	@echo "------------ Penser à: ------------"
	@echo "-----------------------------------"
	@echo "---- changes.txt               ----"
	@echo "---- make test                 ----"
	@echo "---- make VERSION=4.11 version ----"
	@echo "---- cvs commit                ----"
	@echo "---- make VERSION=4.11 rtag    ----"
	@echo "---- make package              ----"
	@echo "---- make submit               ----"
	@echo "-----------------------------------"

version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz doc_ps/$(Package).ps version description.txt \
	   changes.txt $(WWW_dir)
	cvs update $(LOCAL)/Triangulation3
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl teillaud < $(TMP_file)
	rm -f $(TMP_file)

test:
	$(MAKE) -C test/$(Package)

clean:
	$(MAKE) -C doc_tex/basic/$(Package) clean
	$(MAKE) -C examples/$(Package) clean
	$(MAKE) -C demo/$(Package) clean
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.* --exclude=Makefile \
	--exclude=cc_make_ref_pages --exclude=cc_ref_wizard \
	--exclude=TODO \
	--exclude=.pure --exclude=*.pv \
	--exclude=.nfs*\
	--exclude=HTML_TEST \
	version maintainer description.txt changes.txt \
	demo include doc_tex doc_ps src test examples

doc: 
	$(MAKE) -C doc_tex/basic/$(Package) clean
	$(MAKE) -C doc_tex/basic/$(Package) $(Package).ps
	cp doc_tex/basic/$(Package)/$(Package).ps doc_ps/

$(Package).ps: 
	$(MAKE) -C doc_tex/basic/$(Package) $(Package).ps
	cp doc_tex/basic/$(Package)/$(Package).ps doc_ps/$(Package).ps

