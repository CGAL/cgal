WRAPPER=wrapper
PDF_WRAPPER=pdf_wrapper
HTML_TARGET_DIR=html
TOOLS_DOC_SOURCE=/KM/usr3/cgal/Manual_tools/doc/
SUBDIR_NAME=DManual
INSTALL_DIR=/www/pub/www.cgal.org/Members/$(SUBDIR_NAME)

EXT_WRAPPER=wrapper_ext
EXT_PDF_WRAPPER=pdf_wrapper_ext
EXT_HTML_TARGET_DIR=html_ext
EXT_INSTALL_DIR=/www/pub/www.cgal.org/$(SUBDIR_NAME)

internal: clean ps_doc pdf_doc html_doc

external: clean ext_ps_doc ext_pdf_doc ext_html_doc

dvi:
	latex $(WRAPPER)
	makeindex $(WRAPPER)
	index_fix $(WRAPPER).ind
	bibtex $(WRAPPER)
	latex $(WRAPPER)

ext_dvi:  
	latex $(EXT_WRAPPER)
	makeindex $(EXT_WRAPPER)
	index_fix $(EXT_WRAPPER).ind
	bibtex $(EXT_WRAPPER)
	latex $(EXT_WRAPPER)

ps_doc:
	latex $(WRAPPER)
	makeindex $(WRAPPER)
	index_fix $(WRAPPER).ind
	bibtex $(WRAPPER)
	latex $(WRAPPER)
	latex $(WRAPPER)
	dvips -o developers_manual.ps $(WRAPPER)

ext_ps_doc: 
	latex $(EXT_WRAPPER)
	makeindex $(EXT_WRAPPER)
	index_fix $(EXT_WRAPPER).ind
	bibtex $(EXT_WRAPPER)
	latex $(EXT_WRAPPER)
	latex $(EXT_WRAPPER)
	dvips -o developers_manual_ext.ps $(EXT_WRAPPER)

pdf_doc: 
	latex $(PDF_WRAPPER)
	makeindex $(PDF_WRAPPER)
	index_fix $(PDF_WRAPPER).ind
	bibtex $(PDF_WRAPPER)
	latex $(PDF_WRAPPER)
	latex $(PDF_WRAPPER)
	dvips -o pdf_developers_manual.ps $(PDF_WRAPPER)
	ps2pdf pdf_developers_manual.ps
	mv pdf_developers_manual.pdf developers_manual.pdf

ext_pdf_doc:
	latex $(EXT_PDF_WRAPPER)
	makeindex $(EXT_PDF_WRAPPER)
	index_fix $(EXT_PDF_WRAPPER).ind
	bibtex $(EXT_PDF_WRAPPER)
	latex $(EXT_PDF_WRAPPER)
	latex $(EXT_PDF_WRAPPER)
	dvips -o pdf_developers_manual_ext.ps $(EXT_PDF_WRAPPER)
	ps2pdf pdf_developers_manual_ext.ps
	mv pdf_developers_manual_ext.pdf developers_manual_ext.pdf

ext_html_doc:
	cc_manual_to_html -o $(EXT_HTML_TARGET_DIR) -title "CGAL Developers' Manual (2.2.4)" -author "<A HREF=""http://www.cgal.org/CGAL"" target=""_top"">The CGAL project</A>" -date "`date +%b` `date +%e`, `date +%Y`" -aux $(EXT_WRAPPER).aux -bbl $(EXT_WRAPPER).bbl $(EXT_WRAPPER).tex
	cp manual_index.html $(EXT_HTML_TARGET_DIR)
	cp title.html $(EXT_HTML_TARGET_DIR)

html_doc:
	cc_manual_to_html -o $(HTML_TARGET_DIR) -title "CGAL Developers' Manual (2.2.4)" -author "<A HREF=""http://www.cgal.org/CGAL"" target=""_top"">The CGAL project</A>" -date "`date +%b` `date +%e`, `date +%Y`" -aux $(WRAPPER).aux -bbl $(WRAPPER).bbl $(WRAPPER).tex
	cp manual_index.html $(HTML_TARGET_DIR)
	cp title.html $(HTML_TARGET_DIR)

html_title:
	cc_manual_to_html -o title_html -title "CGAL Developers' Manual (2.2.4)" -author "<A HREF=""http://www.cgal.org"" target=""_top"">The CGAL project</A>" -date "`date +%b` `date +%e`, `date +%Y`" -show_main title_wrap.tex

install_internal: install_html install_ps install_pdf install_tools_doc
	chmod -R a+r $(INSTALL_DIR)

install_html:
	mkdir -p $(INSTALL_DIR)
	cp -R $(HTML_TARGET_DIR) $(INSTALL_DIR)
	cd ..; gtar zcvf $(INSTALL_DIR)/developers_manual.html.tar.gz Developers_manual/html/*

install_ps:
	mkdir -p $(INSTALL_DIR)
	cp developers_manual.ps $(INSTALL_DIR)
	cd $(INSTALL_DIR); gzip -f developers_manual.ps

install_pdf:
	mkdir -p $(INSTALL_DIR)
	cp developers_manual.pdf $(INSTALL_DIR)

install_tools_doc:
	mkdir -p $(INSTALL_DIR)
	cp $(TOOLS_DOC_SOURCE)/cc_extract.ps.gz $(INSTALL_DIR)/html
	cp $(TOOLS_DOC_SOURCE)/cc_manual.ps.gz $(INSTALL_DIR)/html
	cp $(TOOLS_DOC_SOURCE)/cc_manual_index.ps.gz $(INSTALL_DIR)/html
	cp $(TOOLS_DOC_SOURCE)/latex_to_html.ps.gz $(INSTALL_DIR)/html

install_external: install_ext_html install_ext_ps install_ext_pdf 
	chmod -R a+r $(EXT_INSTALL_DIR)

install_ext_html:
	mkdir -p $(EXT_INSTALL_DIR)/html
	cp -R $(EXT_HTML_TARGET_DIR)/* $(EXT_INSTALL_DIR)/html
	rm -rf $(EXT_INSTALL_DIR)/html/images/CVS
	cd $(EXT_INSTALL_DIR)/..; gtar zcvf $(EXT_INSTALL_DIR)/developers_manual.html.tar.gz $(SUBDIR_NAME)/html/*

install_ext_ps:
	mkdir -p $(EXT_INSTALL_DIR)
	cp developers_manual_ext.ps $(EXT_INSTALL_DIR)/developers_manual.ps
	cd $(EXT_INSTALL_DIR); gzip -f developers_manual.ps

install_ext_pdf:
	mkdir -p $(EXT_INSTALL_DIR)
	cp developers_manual_ext.pdf $(EXT_INSTALL_DIR)/developers_manual.pdf


clean:
	\rm -f *.aux *.bbl *.blg *.dvi *.idx *.ilg *.ind *.ind.unfixed *.log *.toc *.out
