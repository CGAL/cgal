# Makefile for Circulator Package
PACKAGE=Circulator

WWW=/u/www/0/www/htdocs/files/equipes/prisme
WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Circulator
TMP_file=/tmp/submit-mail

MAINTAINER=Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

package: 
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- make test                ----"
	@echo "---- editer changes.txt       ----"
	@echo "---- cvs update; cvs commit ------"
	@echo "---- more version ----------------"
	@echo "---- VERSION=   ; export VERSION -"
	@echo "---- make Doc --------------------"
	@echo "---- make Version-----------------"
	@echo "---- make tag --------------------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version
	echo  "maintainer : $(MAINTAINER)" >> version

tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"` \
	doc_tex examples include test

submit: 
	rm -rf Circulator.tar.gz
	cd test/Circulator ; make clean
	cd examples/Circulator ; make clean
	tar -zcvf Circulator.tar.gz\
	--exclude=CVS --exclude=wrapper.* --exclude=Makefile \
	--exclude=TODO --exclude=run_testsuite --exclude=error.txt \
	--exclude=circulator_compat.h \
	version description.txt long_description.txt changes.txt \
	include doc_tex doc_ps test examples 
	cp -f Circulator.tar.gz *.txt version $(WWW_dir)
	cp -rf doc_ps/* $(WWW_dir)/doc_ps/
	cp -rf doc_html/* $(WWW_dir)/doc_html/

	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Circulator/Circulator.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
	rm -f $(TMP_file)

Doc :
	cvs commit doc_tex
	cd doc_tex/support/Circulator; \
	latex wrapper.tex; bibtex wrapper.tex; \
	latex wrapper.tex; latex wrapper.tex; \
	dvips wrapper -o; gzip wrapper.ps; rm wrapper.dvi
	mv doc_tex/support/Circulator/wrapper.ps.gz doc_ps/Circulator.ps.gz
	cd doc_tex/support/Circulator; \
	cc_manual_to_html -title Circulator -release $(VERSION) \
	-o  ../../../doc_html/support/Circulator main.tex
	echo "done"