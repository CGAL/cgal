#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=C3
Package=C3

.PHONY: Package Test Submit Version Rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

Package: clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make Test                ----"
	@echo "---- make VERSION=X.X Version ----"
	@echo "---- make Commit              ----"
	@echo "---- make VERSION=X.X Rtag    ----"
	@echo "---- make Package             ----"
	@echo "---- make submit              ----"

	@echo "----------------------------------"

symlinks:
	cd include/CGAL/Cartesian ; sh ./create_symlinks_3

Test:
	cd test/Kernel ; ./cgal_test

Version:
	echo "$(VERSION) (" `date +"%e %b %Y"` ")" > version

Commit:
	cvs commit

Rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

Submit: Package
	cp $(Package).tar.gz version description.txt long_description.txt \
	   changes.txt $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mutt -s autohandle cgal-submit@cs.uu.nl <  $(TMP_file)
	rm -f $(TMP_file)

clean:
	rm -f $(Package).tar.gz
	cd include/CGAL/Cartesian ; sh ./remove_symlinks_3
	cd test/Kernel ; make distclean

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO \
	version long_description.txt description.txt changes.txt \
	include
