# Makefile for Triangulation Package

WWW=/u/www/0/www/htdocs/files/equipes/prisme
WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Triangulation_2
TMP_file=/tmp/submit-mail
CGAL_MAKEFILE=/u/alcor/0/prisme_util/CGAL/CGAL-2.1/make/makefile_i686_Linux-2.2.13_g++-2.95.2_LEDA
MAINTAINER=Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

package: 
	@echo "----------------------------------"
	@echo "------------ Faire :-- -----------"
	@echo "----------------------------------"
	@echo "---- cvs update; -----------------"
	@echo "---- make Test                ----"
	@echo "---- editer changes.txt       ----"
	@echo "---- more version ----------------"
	@echo "---- VERSION=   ; export VERSION -"
	@echo "---- cvs commit ------------------"
	@echo "---- make Doc --------------------"
	@echo "---- make Version-----------------"
	@echo "---- make Tag --------------------"
	@echo "---- make Package ----------------"
	@echo "---- make Submit -----------------"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version 
	echo  "maintainer : $(MAINTAINER)" >> version

Tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"` demo doc examples include test

Package: 
	rm -rf Triangulation.tar.gz
	cd test/Triangulation ; make realclean
	cd examples/Triangulation ; make realclean
	cd demo/Triangulation ; make clean
	tar -zcvf Triangulation.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile* \
	--exclude=Makefilevc \
	--exclude=TODO --exclude=run_testsuite --exclude=error.txt \
	version description.txt changes.txt \
	include doc_tex doc_ps test examples demo

Submit:
	cp -f Triangulation.tar.gz *.txt version $(WWW_dir)
	cp -rf doc_ps/* $(WWW_dir)/doc_ps/
	cp -rf doc_html/* $(WWW_dir)/doc_html/

	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Triangulation_2/Triangulation.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
	rm -f $(TMP_file)

Doc:
	cd /u/alcor/0/prisme_util/CGAL/Local/Triangulation
	cvs commit doc
	rm -rf doc_tex/basic/Triangulation/*.*
	rm -rf doc_ps/triangulation.ps.gz
	cp doc/*.tex doc/*.ltex doc/*.ipe doc/*.gif doc_tex/basic/Triangulation
	cp doc/*.gif doc/*.png doc_html/basic/Triangulation
	rm -rf doc_tex/support/Handle/*.*
	cp doc/pointer.tex doc_tex/support/Handle/main.tex
	cd doc;\
	latex wrapper.tex ; bibtex wrapper.tex; latex wrapper.tex; \
	dvips wrapper -o; gzip wrapper.ps; rm wrapper.dvi
	mv doc/wrapper.ps.gz doc_ps/triangulation.ps.gz
	cd doc; cc_manual_to_html -title Triangulations -release $(VERSION) \
	-o ../doc_html/basic/Triangulation main.tex
	echo "done"

