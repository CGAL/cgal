#  Makefile for CGAL's Triangulation3 package.

#CVS_MODULE_NAME=tet
Package=Viewer_3

#.PHONY: doc package test submit version rtag
VERSION=0.5
WWW_base=/u/www-sop/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: $(Package).tar.gz
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
	echo "maintainer: Francois Rebufat <Francois.Rebufat@sophia.inria.fr>" >> version  

#rtag:
#	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz doc_ps/$(Package).ps version description.txt \
	    $(WWW_dir)
#	cvs update $(LOCAL)/Tetrahedralization
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl rebufat < $(TMP_file)
	rm -f $(TMP_file)

test: ./test/Viewer_3/Viewer_test.o
	$(MAKE) -C test/$(Package)

clean:
#	$(MAKE) -C doc_tex/support/$(Package) clean
#	$(MAKE) -C examples/$(Package) clean
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=manual.tex --exclude=Makefile \
        --exclude=Viewer_3_ref \
	--exclude=TODO \
	--exclude=.pure --exclude=*.pv \
	--exclude=.nfs*\
	version description.txt \
	demo include doc_tex doc_ps src 

doc: doc_ps/$(Package).ps
	$(MAKE) -C doc_tex/support/$(Package)
	cp doc_tex/support/$(Package)/$(Package).ps doc_ps/
