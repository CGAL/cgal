# $Id$
#-----------------------------------------------------------

WWW = /u/zenon/5/www/files/equipes/prisme

FW = fw

.PHONY:	tex

doc:
	latex wrapper ; latex wrapper
	dvips -o wrapper.ps wrapper
	cp geomview.tex  ../doc_tex/support/Geomview/main.tex
	gzip wrapper.ps
	mv wrapper.ps.gz ../doc_ps/geomview.ps.gz
	cc_manual_to_html -title Geomview -release 0.9 \
		-cgal_dir $CGAL/include/ -o ../doc_html/support/Geomview \
		geomview.tex


.PHONY:	print

clean:
	- rm -rf *log *aux *dvi *toc *loc *bbl *blg
	
#-----------------------------  RCS  -----------------------
# $Log$
# Revision 1.1  1998/09/30 12:15:00  hbronni
# Initial revision
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
