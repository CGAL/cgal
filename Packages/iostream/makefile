# makefile created by M Yvinec
# on March 99, 9th
#-----------------------------------------------------------

WWW = /u/www/0/www/htdocs/files/equipes/prisme

VERSION=1.9.1
package:
	echo "Avez vous mis a jour le numero de version"
	echo "Avez vous penser a modifier les fichiers (long_)description.txt"
	echo "Avez vous mis a jour la doc"
	rm -rf iostream.tar.gz
	rm -rf version ; echo "$(VERSION) (`date '+%d %h %y'`)" > version
	tar -zcvf iostream.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=makefile \
	--exclude=TODO \
	version description.txt long_description.txt \
	src include doc_tex doc_ps 
	cp -f iostream.tar.gz *description.txt $(WWW)/CGAL/DR:/Members/IO/
	cp -rf doc_ps/* $(WWW)/CGAL/DR:/Members/IO/doc_ps/
	cp -rf doc_html/* $(WWW)/CGAL/DR:/Members/IO/doc_html/
	echo "Penser a mettre a jour le fichier $(WWW)/CGAL/DR:/Members/index.html"
	echo "Penser a envoyer un mail a cgal-submit quand c'est pret"
	echo "Penser a tagger les fichiers soumis avec le numero de version"

doctex:
	cd /u/alcor/0/prisme_util/CGAL/Local/IO/iostream
	cvs commit doc
	rm -f doc_tex/general/*.*
	rm -f doc_ps/iostream.ps.gz
	cp doc/*.tex doc/*.ltex doc/*.ipe doc_tex/general
	cd doc;\
	latex wrapper.tex ; bibtex wrapper.tex; latex wrapper.tex; \
	dvips wrapper -o iostream.ps; gzip iostream.ps; rm wrapper.dvi
	mv doc/iostream.ps.gz doc_ps
	echo "done"