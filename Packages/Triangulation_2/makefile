# Makefile for Triangulation Package

WWW=/u/www/0/www/htdocs/files/equipes/prisme
WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Triangulation_2
TMP_file=/tmp/submit-mail

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

tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"` demo doc examples include test

submit: 
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
	cp doc/*.gif doc_html/basic/Triangulation
	rm -rf doc_tex/support/Handle/*.*
	cp doc/pointer.tex doc_tex/support/Handle/main.tex
	cd doc;\
	latex wrapper.tex ; bibtex wrapper.tex; latex wrapper.tex; \
	dvips wrapper -o; gzip wrapper.ps; rm wrapper.dvi
	mv doc/wrapper.ps.gz doc_ps/triangulation.ps.gz
	cd doc; cc_manual_to_html -title Triangulations -release $(VERSION) \
	-o ../doc_html/basic/Triangulation main.tex
	echo "done"

.PHONY:	print

print: tex
	dvips -o /tmp/constrained_triangulation.ps constrained_triangulation
	psnup -2 /tmp/constrained_triangulation.ps | lpr 
	rm -f /tmp/constrained_triangulation.ps
	dvips -o /tmp/triangulation.ps triangulation
	psnup -2 /tmp/triangulation.ps | lpr 
	rm -f /tmp/triangulation.ps
	dvips -o /tmp/Pointer.ps Pointer
	psnup -2 /tmp/Pointer.ps | lpr 
	rm -f /tmp/Pointer.ps

#-----------------------------  RCS  -----------------------
# $Log$
# Revision 1.2  2000/01/04 08:10:39  yvinec
# added to CVS files : TODO changes.txt description.txt makefile version
#
# Revision 1.1.4.2  1999/07/08 08:40:23  yvinec
# exclude error.txt from submission
#
# Revision 1.1.4.1  1999/06/21 15:21:17  yvinec
# copy makefile on branch b18juin
#
# Revision 1.1.2.11  1999/05/28 07:36:09  yvinec
#
# include .gif files in doc_tex
#
# Revision 1.1.2.10  1999/05/21 17:20:32  yvinec
# OK
#
# Revision 1.1.2.9  1999/05/12 15:34:29  yvinec
# ok
#
# Revision 1.1.2.8  1999/05/06 12:58:41  yvinec
# OK
#
# Revision 1.1.2.7  1999/04/12 08:40:29  yvinec
# OK
#
# Revision 1.1.2.6  1999/04/09 08:20:01  yvinec
# change make file to make packaage submit more automatic
#
# Revision 1.1.2.5  1999/04/08 15:50:53  yvinec
# changed index.html to shtml in makefile
#
# Revision 1.1.2.4  1999/03/31 14:51:09  yvinec
# reivnclude pre and post conditions
#
# Revision 1.1.2.2  1999/03/09 14:04:19  yvinec
# added conversion to html manual
#
# Revision 1.1.2.1  1999/03/09 13:08:41  yvinec
# makefile incudes
# make package
# make doctex
#
# Revision 1.12  1998/02/18 15:11:18  hbronni
# changed version number
#
# Revision 1.11  1998/02/17 13:37:28  hbronni
# changed version number
#
# Revision 1.10  1998/02/17 12:36:06  hbronni
# added doc_html in target package
#
# Revision 1.9  1998/02/12 18:41:45  hbronni
# added actions in tex / print for constrained_triangulation
#
# Revision 1.8  1998/02/12 12:26:45  hbronni
# added fw +Q constrained_triangulation to target code
#
# Revision 1.7  1998/01/30 13:49:10  hbronni
# Changed triangulation version to 1.12
#
# Revision 1.6  1997/12/18 18:06:00  hbronni
# changed to version 1.11
#
# Revision 1.5  1997/12/17 16:53:38  hbronni
# corrected typo in rcs -n, put rcs -N"Version": instead
#
# Revision 1.4  1997/12/17 16:43:11  hbronni
# Added RCS_VERSION so that files keep track of proper changes
# and appropriate version numbers
#
# Revision 1.3  1997/12/17 16:41:14  hbronni
# changed version number to 1.10 for next release
#
# Revision 1.2  1997/12/16 13:50:08  hbronni
# corrected make package so that it updates the package automatically.
# Usage: check VERSION in the makefile, and type make package.
# Don't forget to send a message to cgal-submit afterwards to notify of the change.
#
# Revision 1.1  1997/12/15 14:34:49  hbronni
# Initial revision
#
#-----------------------------  RCS  -----------------------
