all : doc_ps doc_html

#LATEX_CONV_INPUTS=.:../../examples:Triangulation_2:TDS_2:Triangulation_2_ref:TDS2_ref:
#TEXINPUTS=.:/0/prisme_util/latex:/u/abeille/0/prisme/yvinec/tex/inputs:Triangulation_2:TDS_2:../../examples:Triangulation_2_ref:TD2_ref:

doc_ps :
	latex wrapper
	bibtex wrapper
	makeindex wrapper
	index_fix wrapper.ind
	latex wrapper
	dvips wrapper -o

doc_html :
	rm -rf  ../../doc_html/basic/
	mkdir -p ../../doc_html/basic/Triangulation_2
	mkdir -p ../../doc_html/basic/Triangulation_2_ref
	mkdir -p ../../doc_html/basic/TDS_2
	mkdir -p ../../doc_html/basic/TDS_2_ref
	LATEX_CONV_INPUTS=.:Triangulation_2:TDS_2:Triangulation_2_ref:TDS2_ref\
             cc_manual_to_html -o ../../doc_html/basic/Triangulation_2\
	         Triangulation_2/main.tex
	LATEX_CONV_INPUTS=.:Triangulation_2:TDS_2:Triangulation_2_ref:TDS2_ref\
              cc_manual_to_html  -o ../../doc_html/basic/Triangulation_2_ref \
                 Triangulation_2_ref/main.tex 
	LATEX_CONV_INPUTS=.:Triangulation_2:TDS_2:Triangulation_2_ref:TDS2_ref\
               cc_manual_to_html -o ../../doc_html/basic/TDS_2  \
	         TDS_2/main.tex
	LATEX_CONV_INPUTS=.:Triangulation_2:TDS_2:Triangulation_2_ref:TDS2_ref\
              cc_manual_to_html  -o ../../doc_html/basic/TDS_2_ref \
	         TDS_2_ref/main.tex
	cp Triangulation_2/*.gif ../../doc_html/basic/Triangulation_2
	cp Triangulation_2_ref/*.gif ../../doc_html/basic/Triangulation_2_ref
	cp TDS_2/*.gif ../../doc_html/basic/TDS_2
	cp TDS_2_ref/*.gif ../../doc_html/basic/TDS_2_ref

#
#
#
#
