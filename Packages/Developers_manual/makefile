WRAPPER_FILE=wrapper
HTML_TARGET_DIR=html
EXAMPLE_SOURCE=/KM/guests/hert/cgal/Tools/example
SCRIPTS_SOURCE=/KM/guests/hert/cgal/bin
TOOLS_DOC_SOURCE=/KM/guests/hert/cgal/Tools/doc/
INSTALL_DIR=/www/pub/www.cgal.org/Members/DManual

dvi:
	latex $(WRAPPER_FILE)
	makeindex $(WRAPPER_FILE)
	index_fix $(WRAPPER_FILE).ind
	bibtex $(WRAPPER_FILE)
	latex $(WRAPPER_FILE)

ps_doc:
	latexloop -idx -bib $(WRAPPER_FILE)
	dvips -o developers_manual.ps $(WRAPPER_FILE)
	/opt/adobe/AcroDist_2.1/bin/distill developers_manual.ps

html_doc:
	cc_manual_to_html -o $(HTML_TARGET_DIR) -title "CGAL Developers' Manual (2.1)" -author "<A HREF=""http://www.mpi-sb.mpg.de/GALIA"" target=""_top"">The GALIA project</A>" -date "`date +%b` `date +%e`, `date +%Y`" -aux wrapper.aux -bbl wrapper.bbl $(WRAPPER_FILE).tex
	cp manual_index.html $(HTML_TARGET_DIR)
	cp title.html $(HTML_TARGET_DIR)

html_title:
	cc_manual_to_html -o title_html -title "CGAL Developers' Manual (2.1)" -author "<A HREF=""http://www.mpi-sb.mpg.de/GALIA"" target=""_top"">The GALIA project</A>" -date "`date +%b` `date +%e`, `date +%Y`" -show_main title_wrap.tex

install: install_local install_html install_ps install_pdf install_tools_doc
	chmod -R a+r $(INSTALL_DIR)

install_local: install_html_example install_html_scripts

install_html_example:
	mkdir -p $(HTML_TARGET_DIR)/example/doc_tex/basic/Package/Package_ref
	cp -r $(EXAMPLE_SOURCE)/* $(HTML_TARGET_DIR)/example

install_html_scripts:
	mkdir -p $(HTML_TARGET_DIR)/scripts
	cp $(SCRIPTS_SOURCE)/pstogif $(HTML_TARGET_DIR)/scripts
	cp $(SCRIPTS_SOURCE)/eps2gif $(HTML_TARGET_DIR)/scripts
	cp $(SCRIPTS_SOURCE)/ipe2gif $(HTML_TARGET_DIR)/scripts
	cp $(SCRIPTS_SOURCE)/latex2gif $(HTML_TARGET_DIR)/scripts
	cp $(SCRIPTS_SOURCE)/rename_clib_calls $(HTML_TARGET_DIR)/scripts

install_html:
	cp -R $(HTML_TARGET_DIR) $(INSTALL_DIR)
	cd ..; gtar zcvf $(INSTALL_DIR)/developers_manual.html.tar.gz developers_manual/html/*

install_ps:
	cp developers_manual.ps $(INSTALL_DIR)
	cd $(INSTALL_DIR); gzip -f developers_manual.ps

install_pdf:
	cp developers_manual.pdf $(INSTALL_DIR)

install_tools_doc:
	cp $(TOOLS_DOC_SOURCE)/cc_extract.ps.gz $(INSTALL_DIR)/html
	cp $(TOOLS_DOC_SOURCE)/cc_manual.ps.gz $(INSTALL_DIR)/html
	cp $(TOOLS_DOC_SOURCE)/cc_manual_index.ps.gz $(INSTALL_DIR)/html
	cp $(TOOLS_DOC_SOURCE)/latex_to_html.ps.gz $(INSTALL_DIR)/html

clean:
	latexclean
