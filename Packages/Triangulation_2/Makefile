
# Makefile for Triangulation Package

PACKAGE=Triangulation

WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Triangulation_2
TMP_file=/tmp/submit-mail
CGAL_MAKEFILE=/u/alcor/0/prisme_util/CGAL/CGAL-I/make/makefile_i686_Linux-2.2.14_g++-2.95.2_LEDA
TEXINPUTS=.:../../../examples:/0/prisme_util/latex:/u/procyon/0/prisme/yvinec/tex/inputs:
LATEX_CONV_INPUTS=.:../../../examples:/u/procyon/0/prisme/yvinec/Cgal/Triangulation/examples:

package: 
	@echo "----------------------------------"
	@echo "------------ Faire :-- -----------"
	@echo "----------------------------------"
	@echo "---- cvs update; -----------------"
	@echo "---- make Test ; ------------ ----"
	@echo "---- make VERSION=4.84 Doc -------"
	@echo "---- editer changes.txt       ----"
	@echo "---- make VERSION=4.84 Version ---"
	@echo "---- cvs commit ------------------"
	@echo "---- make VERSION=4.84 Tag -------"
	@echo "---- make VERSION=4.84 Submit ----"
	@echo "---------------------------------"


Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version 

Tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"` \
	demo doc_tex examples include test


Test: 
	cd test/Triangulation2_Delaunay; mv makefile makefile.submitted
	cd test/Triangulation2_constrained; mv makefile makefile.submitted
	cd test; ./run_testsuite;
	cd test/Triangulation; make realclean
	cd test/Triangulation2_Delaunay; make clean
	cd test/Triangulation2_constrained; make clean
	cd test/Triangulation2_newd; make -f Makefile realclean;
	cd test/Triangulation2_Delaunay; mv makefile.submitted makefile
	cd test/Triangulation2_constrained; mv makefile.submitted makefile


clean:
	rm -f $(PACKAGE).tar.gz
	cd test/Triangulation ; make -f Makefile realclean;
	cd test/Triangulation2_Delaunay; make -f Makefile realclean
	cd test/Triangulation2_constrained; make -f Makefile realclean
	cd test/Triangulation2_newd; make -f Makefile realclean;
	cd examples/Triangulation ; make realclean
	cd demo/Triangulation ; make realclean


$(PACKAGE).tar.gz:
	tar -zcvf $(PACKAGE).tar.gz \
	--exclude=CVS --exclude=wrapper.* --exclude=Makefile* \
	--exclude=Makefilevc \
	--exclude=TODO --exclude=run_testsuite --exclude=error.txt \
	maintainer version description.txt changes.txt \
	include doc_tex doc_ps test examples demo

Submit: clean $(PACKAGE).tar.gz
	cp -f $(PACKAGE).tar.gz *.txt version $(WWW_dir)
	cp -rf doc_ps/* $(WWW_dir)/doc_ps/
	cp -rf doc_html/* $(WWW_dir)/doc_html/

	cgal_submit $(VERSION)

#	echo > $(TMP_file)
#	echo "submission::" >> $(TMP_file)
#	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
#	echo "/Triangulation_2/Triangulation.tar.gz" >> $(TMP_file)
#	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
#	rm -f $(TMP_file)

Doc : Doc_ps Doc_html


Doc_ps:
	rm -f doc_ps/Triangulation.ps.gz;
	cd doc_tex/basic/Triangulation;\
	latex wrapper.tex; bibtex wrapper; \
	latex wrapper.tex; latex wrapper.tex; \
	dvips -o Triangulation.ps wrapper.dvi; \
	gzip Triangulation.ps;\
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl *.toc ;
	mv doc_tex/basic/Triangulation/Triangulation.ps.gz doc_ps ; 
	rm -f doc_ps/pointer.ps.gz;
	cd doc_tex/support/Handle;\
	latex wrapper.tex; bibtex wrapper; \
	latex wrapper.tex; latex wrapper.tex; \
	dvips -o pointer.ps wrapper.dvi; \
	gzip pointer.ps; \
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl *.toc
	mv doc_tex/support/Handle/pointer.ps.gz doc_ps ; \


Doc_html :
	cd doc_tex/basic/Triangulation; \
	cc_manual_to_html -title Triangulations \
	-release $(VERSION) \
	-o ../../../doc_html/basic/Triangulation main.tex
	cd doc_tex/support/Handle; \
	cc_manual_to_html -title Handles \
	-release $(VERSION) \
	-o ../../../doc_html/support/Handle main.tex
	echo "done"



