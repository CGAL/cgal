#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=Cartesian_basic
Package=Cartesian_basic

.PHONY: package test submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make Test                ----"
	@echo "---- make VERSION=3.1 version ----"
	@echo "---- make Commit               ----"
	@echo "---- make VERSION=3.1 Rtag    ----"
	@echo "---- make Package             ----"
	@echo "---- make Submit              ----"
	@echo "----------------------------------"

version:
	echo "$(VERSION) (" `date +"%e %b %Y"` ")" > version

commit:
	cvs commit

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz version description.txt changes.txt \
	   $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl < $(TMP_file)
	rm -f $(TMP_file)

clean:
	rm -f $(Package).tar.gz

$(Package).tar.gz:
	./check_headers.pl -u include/CGAL/*.[Ch]
	./check_headers.pl -u include/CGAL/*/*.[Ch]
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO \
	version long_description.txt description.txt changes.txt \
	include
	./check_headers.pl -u include/CGAL/*.[Ch]
	./check_headers.pl -u include/CGAL/*/*.[Ch]
