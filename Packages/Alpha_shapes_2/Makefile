#  Makefile for CGAL's Alpha_shapes_2 package. [Frank DA]

CVS_MODULE_NAME=alpha
Package=Alpha_shapes_2

.PHONY: doc package filters test example submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: mode doc clean $(Package).tar.gz
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

mode:
	chmod a+x include include/CGAL doc_tex doc_tex/basic \
	doc_tex/basic/$(Package) doc_ps demo demo/$(Package) \
	demo/$(Package)/Data examples examples/$(Package) \
	examples/$(Package)/data test test/$(Package) test/$(Package)/data
	chmod -R a+r *	

version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version
	echo "Maintainer : Frank Da <Frank.Da@sophia.inria.fr>" >> version

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz doc_ps/$(Package).ps version description.txt changes.txt \
	   $(WWW_dir)
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	chmod a+r $(WWW_dir)/$(Package).ps.gz
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl frankda < $(TMP_file)
	rm -f $(TMP_file)

test:
	$(MAKE) -C test/$(Package) test_alpha test_weight

example:
	$(MAKE) -C examples/$(Package) example_alpha example_weight

clean:
	$(MAKE) -C doc_tex/basic/$(Package) clean
	$(MAKE) -C test/$(Package) clean
	$(MAKE) -C examples/$(Package) clean
	$(MAKE) -C demo/$(Package) -f Makefile clean	
	rm -f $(Package).tar.gz

# --exclude=tstlink is here because I don't want to provide a makefile for
# the test-suite.
# --exclude=tst6.C because it's supposed to warn at compile time, but it
# actually breaks at run time.

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS \
	--exclude=wrapper.tex --exclude=weight.tex \
	--exclude=Weighted_alpha_shape_2.h \
	--exclude=Makefile* --exclude=*.o \
        --exclude=*.obj --exclude=*.exe --exclude=*~ \
	--exclude=*5000* --exclude=Attali_7 \
	--exclude=Melkemi_1 --exclude=allerton \
	--exclude=baguette* --exclude=test2 --exclude=test.franki \
	version description.txt changes.txt \
	include doc_tex doc_ps test examples demo

doc: doc_ps/$(Package).ps

doc_ps/$(Package).ps: doc_tex/basic/$(Package)/alpha.tex
	$(MAKE) -C doc_tex/basic/$(Package) all
	cp doc_tex/basic/$(Package)/$(Package).ps doc_ps/