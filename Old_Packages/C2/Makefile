#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=C2
Package=C2

.PHONY: package test submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make test                ----"
	@echo "---- make VERSION=3.1 version ----"
	@echo "---- make commit               ----"
	@echo "---- make VERSION=3.1 rtag    ----"
	@echo "---- make package             ----"
	@echo "---- make submit              ----"
	@echo "----------------------------------"

test:
	cd test/Kernel; ./cgal_test

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
	mutt -s autohandle cgal-submit@cs.uu.nl hbronni < $(TMP_file)
	rm -f $(TMP_file)

clean:
	rm -f $(Package).tar.gz
	cd include/CGAL/Cartesian ; sh ./remove_symlinks
	cd test/Kernel ; make clean

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO \
	version long_description.txt description.txt changes.txt \
	include
