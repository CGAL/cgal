#located in doc_tex;

VERSION_NUM=2.4
VERSION=CGAL-$(VERSION_NUM)

# location of the test manul results
TEST_DESTINATION=/www/pub/www.cgal.org/Members/Manual_test
# subdirectory containing the test manual pieces and log files
TEST_SUBDIR=test_manual

# location of the final manual
#DESTINATION=/www/pub/www.cgal.org/Members/test_release/Manual
DESTINATION=/www/pub/www.cgal.org/Manual
# location of the doc_ps subdirectory
DOC_PS_DESTINATION=$(DESTINATION)
# location of the doc_pdf subdirectory
DOC_PDF_DESTINATION=$(DESTINATION)
# location of the doc_html subdirectory
DOC_HTML_DESTINATION=$(DESTINATION)
# location of the tar files for downloading
DOC_DOWNLOAD_DIR=/ftp/pub/outgoing/CGAL
#DOC_DOWNLOAD_DIR=/www/pub/www.cgal.org/Members/test_release/Manual

# location of the frame set and the extra index files of the manual
HTML_EXTRA_SOURCE=/KM/usr3/cgal/web_pages/Manual
# location of the frame set and the extra index files of the manual
INCLUDE_SOURCE=/KM/projects/CGAL/Work_in_Progress/releases/CGAL-2.4/

all:
	driver
	\rm -f basic/*/*.aux 
	\rm -f support/*/*.aux 
	to_pdf

install_all: install pack dest_zip_ps
	
install: title_html doc_html doc_ps doc_pdf install_html_extras
	mkdir -p $(DESTINATION)
	cp -r $(INCLUDE_SOURCE)/include $(DESTINATION)
	chmod -R a+rX $(DESTINATION)

test: 
	test_driver
	create_test_page $(VERSION) $(TEST_SUBDIR)

test_page:
	create_test_page $(VERSION) $(TEST_SUBDIR)

test_basic: 
	test_driver basic;

test_support: 
	test_driver support;

test_general: 
	test_driver general;

test_kernel: 
	test_driver kernel;

test_kernel_d: 
	test_driver kernel_d;

test_installation: 
	test_driver installation;

install_test: title_html install_test_page install_log_files \
              install_test_ps install_test_html

pack: pack_doc_ps_a4 pack_doc_ps_letter pack_doc_pdf pack_doc_html

pack_doc_ps_a4:
	cd $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4; \
           mkdir doc_ps; cp *.ps* doc_ps; \
           cp ../index.html doc_ps; \
           tar -cvf $(VERSION)-doc_ps.tar doc_ps ; \
           gzip -f $(VERSION)-doc_ps.tar; \
           gtar -cvf $(VERSION)-doc_ps.gtar doc_ps ; \
           gzip -f $(VERSION)-doc_ps.gtar; \
           zip -r $(VERSION)-doc_ps.zip doc_ps ; \
           mv $(VERSION)-doc_ps.*.gz $(DOC_DOWNLOAD_DIR); \
           mv $(VERSION)-doc_ps.zip $(DOC_DOWNLOAD_DIR); \
           \rm -r doc_ps 
	chmod a+r $(DOC_DOWNLOAD_DIR)/$(VERSION)-doc_ps*;

pack_doc_ps_letter:
	cd $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter;\
           mkdir doc_ps; cp *.ps* doc_ps; \
           cp ../index.html doc_ps; \
           tar -cvf $(VERSION)-doc_ps_letter.tar doc_ps ; \
           gzip -f $(VERSION)-doc_ps_letter.tar; \
           gtar -cvf $(VERSION)-doc_ps_letter.gtar doc_ps ; \
           gzip -f $(VERSION)-doc_ps_letter.gtar; \
           zip -r $(VERSION)-doc_ps_letter.zip doc_ps ; \
           mv $(VERSION)-doc_ps_letter.*.gz $(DOC_DOWNLOAD_DIR); \
           mv $(VERSION)-doc_ps_letter.zip $(DOC_DOWNLOAD_DIR); \
           \rm -r doc_ps
	chmod a+r $(DOC_DOWNLOAD_DIR)/$(VERSION)-doc_ps*;

pack_doc_pdf:
	cd $(DOC_PDF_DESTINATION)/doc_pdf;\
           mkdir doc_pdf; cp *.pdf doc_pdf; \
           cp index.html doc_pdf; \
           tar -cvf $(VERSION)-doc_pdf.tar doc_pdf ; \
           gzip -f $(VERSION)-doc_pdf.tar; \
           gtar -cvf $(VERSION)-doc_pdf.gtar doc_pdf ; \
           gzip -f $(VERSION)-doc_pdf.gtar; \
           zip -r $(VERSION)-doc_pdf.zip doc_pdf ; \
           mv $(VERSION)-doc_pdf.*.gz $(DOC_DOWNLOAD_DIR); \
           mv $(VERSION)-doc_pdf.zip $(DOC_DOWNLOAD_DIR); \
           \rm -r doc_pdf
	chmod a+r $(DOC_DOWNLOAD_DIR)/$(VERSION)-doc_pdf*;

pack_doc_html:
	cd $(DOC_HTML_DESTINATION); \
          tar -cvf $(VERSION)-doc_html.tar doc_html ; \
          gzip -f $(VERSION)-doc_html.tar; \
          gtar -cvf $(VERSION)-doc_html.gtar doc_html ; \
          gzip -f $(VERSION)-doc_html.gtar; \
          zip -r $(VERSION)-doc_html.zip doc_html ; \
          mv $(VERSION)-doc_html.* $(DOC_DOWNLOAD_DIR); \
	chmod a+r $(DOC_DOWNLOAD_DIR)/$(VERSION)-doc_html.*

installation_update:
	driver installation
	to_pdf installation; 
	cd installation; \
          gzip -f installation.ps; \
          cp installation.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/; \
          gzip -f installation_letter.ps; \
          cp installation_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/installation.ps.gz ;\
          cp installation.pdf $(DOC_PDF_DESTINATION)/doc_pdf/
	cd installation/html; \
          cp -R * $(DOC_HTML_DESTINATION)/doc_html/installation/.
	chmod a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/installation.ps.gz
	chmod a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/installation.ps.gz
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/installation

general_update:
	driver general;
	to_pdf general;
	cd general; \
           gzip -f general_intro.ps;  \
           cp general_intro.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/. ;\
           gzip -f general_intro_letter.ps; \
           cp general_intro_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/general_intro.ps.gz ;\
           cp general_intro.pdf $(DOC_PDF_DESTINATION)/doc_pdf ;
	cd general/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/general_intro/.
	chmod -R a+r $(DOC_PS_DESTINATION)/doc_ps
	chmod a+r $(DOC_PDF_DESTINATION)/doc_pdf/general_intro.pdf
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/general_intro

kernel_update:
	driver kernel;
	to_pdf kernel;
	cd kernel; \
           gzip -f kernel.ps; \
           cp kernel.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/. ;\
           gzip -f kernel_letter.ps; \
           cp kernel_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/kernel.ps.gz; \
           cp kernel.pdf $(DOC_PDF_DESTINATION)/doc_pdf 
	cd kernel/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/kernel/.
	chmod -R a+r $(DOC_PS_DESTINATION)/doc_ps/
	chmod a+r $(DOC_PDF_DESTINATION)/doc_pdf/kernel.pdf
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/kernel

kernel_d_update:
	driver kernel_d;
	to_pdf kernel_d;
	cd kernel_d; \
           gzip -f kernel_d.ps; \
           cp kernel_d.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/. ;\
           gzip -f kernel_d_letter.ps; \
           cp kernel_d_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/kernel_d.ps.gz; \
           cp kernel_d.pdf $(DOC_PDF_DESTINATION)/doc_pdf 
	cd kernel_d/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/kernel_d/.
	chmod -R a+r $(DOC_PS_DESTINATION)/doc_ps/
	chmod a+r $(DOC_PDF_DESTINATION)/doc_pdf/kernel_d.pdf
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/kernel_d

deprecated_basic_update:
	driver deprecated_basic;
	to_pdf deprecated_basic;
	cd deprecated_basic; \
           gzip -f basic_lib.ps; \
           cp basic.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/deprecated_basic_lib.ps.gz. ;\
           gzip -f basic_lib_letter.ps; \
           cp basic_lib_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/deprecated_basic_lib_letter.ps.gz; \
           cp basic.pdf $(DOC_PDF_DESTINATION)/doc_pdf/deprecated_basic_lib.pdf
	cd deprecated_basic/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/deprecated_basic_lib/.
	chmod -R a+r $(DOC_PS_DESTINATION)/doc_ps/
	chmod a+r $(DOC_PDF_DESTINATION)/doc_pdf/deprecated_basic_lib.pdf
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/deprecated_basic_lib

basic_update:
	driver basic;
	to_pdf basic; 
	cd basic; \
           gzip -f basic_lib.ps; \
           cp basic_lib.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/. ; \
           gzip -f basic_lib_letter.ps; \
           cp basic_lib_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/basic_lib.ps.gz ; \
	   cp basic_lib.pdf $(DOC_PDF_DESTINATION)/doc_pdf
	cd basic/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/basic_lib/.
	chmod a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/basic_lib.ps.gz
	chmod a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/basic_lib.ps.gz
	chmod a+r $(DOC_PS_DESTINATION)/doc_pdf/basic_lib.pdf
	chmod -R a+rX  $(DOC_HTML_DESTINATION)/doc_html/basic_lib

support_update:
	driver support;
	to_pdf support; 
	cd support; \
           gzip -f support_lib.ps; \
           cp support_lib.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/. ; \
	   gzip -f support_lib_letter.ps; \
           cp support_lib_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/support_lib.ps.gz ; \
	   cp support_lib.pdf $(DOC_PDF_DESTINATION)/doc_pdf/ 
	cd support/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/support_lib/.
	chmod a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/support_lib.ps.gz
	chmod a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/support_lib.ps.gz
	chmod a+r $(DOC_PDF_DESTINATION)/doc_pdf/support_lib.pdf
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/support_lib

use_of_stl_update:
	driver use_of_stl;
	to_pdf use_of_stl;
	cd use_of_stl; \
           gzip -f use_of_stl.ps; \
           cp use_of_stl.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4 ; \
	   gzip -f use_of_stl_letter.ps; \
           cp use_of_stl_letter.ps.gz $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letteruse_of_stl.ps.gz ; \
	   cp use_of_stl.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. 
	cd use_of_stl/html; \
           cp -R * $(DOC_HTML_DESTINATION)/doc_html/use_of_stl/.
	chmod -$ a+r $(DOC_PS_DESTINATION)/doc_ps
	chmod a+r $(DOC_PDF_DESTINATION)/doc_pdf/use_of_stl.pdf
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html/use_of_stl

doc_ps: doc_ps_a4 doc_ps_letter

doc_ps_a4:
	mkdir -p $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps 
	cd installation; cp installation.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;
	cd kernel; cp kernel.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;
	cd kernel_d; cp kernel_d.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;
	cd basic; cp basic_lib.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;
	cd deprecated_basic; cp basic_lib.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/deprecated_basic_lib.ps ;
	cd support; cp support_lib.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;
	cd use_of_stl; cp use_of_stl.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;
	mv $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4
	rmdir $(DOC_PS_DESTINATION)/doc_ps/a4
	chmod -R a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4;

#	cd general; cp general_intro.ps $(DOC_PS_DESTINATION)/doc_ps/a4/doc_ps/. ;

doc_ps_letter:
	mkdir -p $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps
	cd installation; cp installation_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/installation.ps ;
	cd kernel; cp kernel_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/kernel.ps ;
	cd kernel_d; cp kernel_d_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/kernel_d.ps ;
	cd basic; cp basic_lib_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/basic_lib.ps ;
	cd deprecated_basic; cp basic_lib_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/deprecated_basic_lib.ps ;
	cd support; cp support_lib_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/support_lib.ps ;
	cd use_of_stl; cp use_of_stl_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/use_of_stl.ps ;
	mv $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter
	rmdir $(DOC_PS_DESTINATION)/doc_ps/letter
	chmod -R a+r $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter;

#	cd general; cp general_intro_letter.ps $(DOC_PS_DESTINATION)/doc_ps/letter/doc_ps/general_intro.ps ;

doc_pdf:
	mkdir -p $(DOC_PDF_DESTINATION)/doc_pdf
	cd installation; cp installation.pdf $(DOC_PDF_DESTINATION)/doc_pdf/ ;
	cd kernel; cp kernel.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. ;
	cd kernel_d; cp kernel_d.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. ;
	cd basic; cp basic_lib.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. ;
	cd deprecated_basic; cp basic_lib.pdf $(DOC_PDF_DESTINATION)/doc_pdf/deprecated_basic_lib.pdf ;
	cd support; cp support_lib.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. ;
	cd use_of_stl; cp use_of_stl.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. ;
	chmod -R a+r $(DOC_PDF_DESTINATION)/doc_pdf;
#	cd general; cp general_intro.pdf $(DOC_PDF_DESTINATION)/doc_pdf/. ;


zip_ps:
	gzip -f general/general_intro*.ps
	gzip -f kernel/kernel*.ps
	gzip -f kernel_d/kernel_d*.ps
	gzip -f basic/basic_lib*.ps
	gzip -f deprecated_basic/basic_lib*.ps
	gzip -f support/support_lib*.ps
	gzip -f installation/installation*.ps
	gzip -f use_of_stl/use_of_stl*.ps
#	gzip -f general/general_intro*.ps

unzip_ps:
	gunzip kernel/kernel*.ps.gz
	gunzip kernel_d/kernel_d*.ps.gz
	gunzip basic/basic_lib*.ps.gz
	gunzip deprecated_basic/basic_lib*.ps.gz
	gunzip support/support_lib*.ps.gz
	gunzip installation/installation*.ps.gz
	gunzip use_of_stl/use_of_stl*.ps.gz
#	gunzip general/general*.ps.gz

dest_unzip_ps:
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/kernel.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/kernel.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/kernel_d.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/kernel_d.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/basic_lib.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/basic_lib.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/deprecated_basic_lib.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/deprecated_basic_lib.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/support_lib.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/support_lib.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/installation*.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/installation*.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/use_of_stl*.ps.gz
	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/use_of_stl*.ps.gz
#	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/general_intro.ps.gz
#	gunzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/general_intro.ps.gz

dest_zip_ps:
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/kernel.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/kernel.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/kernel_d.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/kernel_d.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/basic_lib.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/basic_lib.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/deprecated_basic_lib.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/deprecated_basic_lib.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/support_lib.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/support_lib.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/installation*.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/installation*.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/use_of_stl*.ps
	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/use_of_stl*.ps
#	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_a4/general_intro.ps
#	gzip -f $(DOC_PS_DESTINATION)/doc_ps/doc_ps_letter/general_intro.ps

install_ps:
	cd installation; cp installation.ps $(DESTINATION)/doc_ps/doc_ps_a4/ ;
	cd installation; cp installation_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/installation.ps ;
	cd kernel; cp kernel.ps $(DESTINATION)/doc_ps/doc_ps_a4/. ;
	cd kernel; cp kernel_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/kernel.ps ;
	cd kernel_d; cp kernel_d.ps $(DESTINATION)/doc_ps/doc_ps_a4/. ;
	cd kernel_d; cp kernel_d_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/kernel_d.ps ;
	cd basic; cp basic_lib.ps $(DESTINATION)/doc_ps/doc_ps_a4/. ;
	cd basic; cp basic_lib_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/basic_lib.ps ;
	cd deprecated_basic; cp basic_lib.ps $(DESTINATION)/doc_ps/doc_ps_a4/deprecated_basic_lib.ps. ;
	cd deprecated_basic; cp basic_lib_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/deprecated_basic_lib.ps ;
	cd support; cp support_lib.ps $(DESTINATION)/doc_ps/doc_ps_a4/ ;
	cd support; cp support_lib_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/support_lib.ps ;
	cd use_of_stl; cp use_of_stl.ps $(DESTINATION)/doc_ps/doc_ps_a4/. ;
	cd use_of_stl; cp use_of_stl_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/use_of_stl.ps ;
#	cd general; cp general_intro.ps $(DESTINATION)/doc_ps/doc_ps_a4 ;
#	cd general; cp general_intro_letter.ps $(DESTINATION)/doc_ps/doc_ps_letter/general_intro.ps ;


title_html:
	cd kernel; cp title.html html/. ;
	cd kernel_d; cp title.html html/. ;
	cd basic; cp Title/title.html html/. ;
	cd deprecated_basic; cp Title/title.html html/. ;
	cd support; cp Title/title.html html/. ;
	cd installation; cp title.html html/. ;
	cd use_of_stl; cp title.html html/. ;

#	cd general; cp title.html html/. ;

doc_html: 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/installation; 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/kernel; 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/kernel_d; 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/basic_lib; 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/deprecated_basic_lib; 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/support_lib; 
	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/use_of_stl; 
	cd kernel/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/kernel/.
	cd kernel_d/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/kernel_d/.
	cd basic/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/basic_lib/.
	cd deprecated_basic/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/deprecated_basic_lib/.
	cd support/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/support_lib/.
	cd installation/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/installation/.
	cd use_of_stl/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/use_of_stl/.
	chmod -R a+rX $(DOC_HTML_DESTINATION)/doc_html

#	mkdir -p $(DOC_HTML_DESTINATION)/doc_html/general_intro; 
#	cd general/html; cp -R * $(DOC_HTML_DESTINATION)/doc_html/general_intro/.

install_html:
	cd installation/html; cp -R * $(DESTINATION)/doc_html/installation/. ;
	cd kernel/html; cp -R * $(DESTINATION)/doc_html/kernel/. ;
	cd kernel_d/html; cp -R * $(DESTINATION)/doc_html/kernel_d/. ;
	cd basic/html; cp -R * $(DESTINATION)/doc_html/basic_lib/. ;
	cd deprecated_basic/html; cp -R * $(DESTINATION)/doc_html/deprecated_basic_lib/. ;
	cd support/html; cp -R * $(DESTINATION)/doc_html/support_lib/. ;
	cd use_of_stl/html; cp -R * $(DESTINATION)/doc_html/use_of_stl/. ;
	chmod -R a+r $(DESTINATION)/doc_html ;
#	cd general/html; cp -R * $(DESTINATION)/doc_html/general_intro/. ;


install_test_page:
	mkdir -p $(TEST_DESTINATION)
	cp $(VERSION).manual_test.html $(TEST_DESTINATION)
	\cp $(VERSION).manual_test.html $(TEST_DESTINATION)/index.html
	chmod a+r $(TEST_DESTINATION)/$(VERSION).manual_test.html
	chmod a+r $(TEST_DESTINATION)/index.html

install_log_files:
	mv $(TEST_DESTINATION)/$(TEST_SUBDIR) $(TEST_DESTINATION)/$(TEST_SUBDIR).bck
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)
	copy_log_files $(TEST_DESTINATION)/$(TEST_SUBDIR)
	chmod -R a+rX $(TEST_DESTINATION)/$(TEST_SUBDIR)
	\rm -r $(TEST_DESTINATION)/$(TEST_SUBDIR).bck

install_test_html: 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/installation; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/kernel; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/kernel_d; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/basic_lib; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/deprecated_basic_lib; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/support_lib; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/use_of_stl; 
	cd installation/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/installation/. ;
	cd kernel/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/kernel/. ;
	cd kernel_d/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/kernel_d/. ;
	cd basic/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/basic_lib/. ;
	cd deprecated_basic/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/deprecated_basic_lib/. ;
	cd support/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/support_lib/. ;
	cd use_of_stl/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/use_of_stl/. ;
	cd $(TEST_DESTINATION)/$(TEST_SUBDIR)/html; gif_grep ;
	chmod -R a+rX $(TEST_DESTINATION)/$(TEST_SUBDIR)/html ;
#	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/general_intro; 
#	cd general/html; cp -R * $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/general_intro/. ;


install_test_ps:
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/installation; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/kernel; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/kernel_d; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/basic_lib; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/deprecated_basic_lib; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/support_lib; 
	mkdir -p $(TEST_DESTINATION)/$(TEST_SUBDIR)/html/use_of_stl; 
	cd installation; gzip -f installation*.ps; cp installation*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/installation/.;
	cd kernel; gzip -f kernel*.ps; cp kernel*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/kernel/. ;
	cd kernel_d; gzip -f kernel_d*.ps; cp kernel_d*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/kernel_d/. ;
	cd basic; gzip -f basic_lib*.ps; cp basic_lib*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/basic/.;
	cd deprecated_basic; gzip -f basic_lib*.ps; cp basic_lib*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/deprecated_basic/.;
	cd support; gzip -f support_lib*.ps; cp support_lib*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/support/. ;
	cd use_of_stl; gzip -f use_of_stl*.ps; cp use_of_stl*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/use_of_stl/. ;
	chmod a+r $(TEST_DESTINATION)/$(TEST_SUBDIR)/*/*.ps.gz
#	cd general; gzip -f general_intro*.ps; cp general_intro*.ps.gz $(TEST_DESTINATION)/$(TEST_SUBDIR)/general/.;

install_html_extras:
	mkdir -p $(DOC_HTML_DESTINATION)
	cp -r $(HTML_EXTRA_SOURCE)/index.html $(DOC_HTML_DESTINATION)
	cp -r $(HTML_EXTRA_SOURCE)/doc_html $(DOC_HTML_DESTINATION)
	rm -rf $(DOC_HTML_DESTINATION)/doc_html/CVS
	cp -r $(HTML_EXTRA_SOURCE)/doc_pdf $(DOC_HTML_DESTINATION)
	rm -rf $(DOC_HTML_DESTINATION)/doc_pdf/CVS
	cp -r $(HTML_EXTRA_SOURCE)/doc_ps $(DOC_HTML_DESTINATION)
	rm -rf $(DOC_HTML_DESTINATION)/doc_ps/CVS

clean:
	rm -f *.aux docdirs_failed* docdirs_passed*;
	cd installation; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl *.ind ionly
	cd general; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl *.ind ionly
	cd kernel; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl *.ind ionly
	cd kernel_d; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl *.ind ionly
	cd basic; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl  *.ind ionly docdirs_failed* docdirs_passed*
	cd deprecated_basic; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl  *.ind ionly docdirs_failed* docdirs_passed*
	cd support; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl *.ind ionly docdirs_failed* docdirs_passed*
	cd use_of_stl; rm -rf *.bck *.aux *.unfixed *.blg *.dvi *.idx *.ilg *.log *.toc *.bbl *.ind ionly
	rm -f installation/installation*.ps.gz
	rm -f general/general_intro*.ps.gz
	rm -f kernel/kernel*.ps.gz
	rm -f kernel_d/kernel_d*.ps.gz
	rm -f basic/basic_lib*.ps.gz
	rm -f deprecated_basic/basic_lib*.ps.gz
	rm -f support/support_lib*.ps.gz
	rm -f use_of_stl/use_of_stl*.ps.gz
