#  Makefile for CGAL's Alpha_shapes_2 package. [Frank DA]

CVS_MODULE_NAME=alpha
Package=Alpha_shapes_2

.PHONY: doc package filters test example submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: clean $(Package).tar.gz
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
	cp $(Package).tar.gz doc_ps/$(Package).ps version description.txt changes.txt \
	   $(WWW_dir)
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl frankda < $(TMP_file)
	rm -f $(TMP_file)

test:
	$(MAKE) -C test/$(Package)

clean:
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz
	rm -f src/.pure src/*.o

# --exclude=tstlink is here because I don't want to provide a makefile for
# the test-suite.
# --exclude=tst6.C because it's supposed to warn at compile time, but it
# actually breaks at run time.

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=*.tex --exclude=Makefile \
	--exclude=TODO --exclude=*5000* --exclude=tstlink2.C \
	--exclude=tst6.C \
	version description.txt changes.txt \
	include doc_tex doc_ps test examples demo


