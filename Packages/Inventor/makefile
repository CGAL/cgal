WWW = /u/www/0/www/htdocs/files/equipes/prisme

VERSION=2.0
package:
	echo "Avez vous mis a jour le numero de version"
	echo "Avez vous penser a modifier les fichiers (long_)description.txt"
	echo "Avez vous mis a jour la doc"
	rm -rf Inventor.tar.gz
	rm -rf version ; echo "$(VERSION) (`date '+%d %h %y'`)" > version
	tar -zcvf Inventor.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=makefile \
	--exclude=TODO --exclude=run_testsuite\
	version description.txt long_description.txt \
	include 
	cp -f Inventor.tar.gz *description.txt $(WWW)/CGAL/DR:/Members/Inventor/
	echo "Penser a mettre a jour le fichier $WWW/CGAL/DR:/Members/index.shtml"
	echo "Penser a envoyer un mail a cgal-submit quand c'est pret"
	echo "Penser a tagger les fichiers soumis avec le numero de version"
	echo "cvs tag v2_x  version description.txt long_description.txt include"
	echo "done"




