#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=tet
Package=Triangulation3

.PHONY: doc package test submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
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
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=Oriented_side_power_test_3.h \
	--exclude=Regular_triangulation_3.h \
	--exclude=Regular_triangulation_geom_traits_3.h \
	--exclude=Weighted_point_3.h \
	--exclude=TODO \
	--exclude=.nfs*\
	version description.txt changes.txt \
	demo include doc_tex doc_ps test examples

doc: doc_ps/$(Package).ps
	$(MAKE) -C doc_tex/basic/$(Package)
	cp doc_tex/basic/$(Package)/$(Package).ps doc_ps/
