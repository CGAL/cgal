all : doc_ps doc_html

LOCALINPUTS=.:../../examples:Triangulation_3:TriangulationDS_3:Triangulation_3_ref:TriangulationDS_3_ref:

doc_ps :
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} latex wrapper
	bibtex wrapper
	makeindex wrapper
	index_fix wrapper.ind
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} latex wrapper
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} dvips wrapper -o

doc_html :
	rm -rf  ../../doc_html/basic/
	mkdir -p ../../doc_html/basic/Triangulation_3
	mkdir -p ../../doc_html/basic/Triangulation_3_ref
	mkdir -p ../../doc_html/basic/TriangulationDS_3
	mkdir -p ../../doc_html/basic/TriangulationDS_3_ref
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
             cc_manual_to_html -o ../../doc_html/basic/Triangulation_3\
	         Triangulation_3/main.tex
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
              cc_manual_to_html  -o ../../doc_html/basic/Triangulation_3_ref \
                 Triangulation_3_ref/main.tex 
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
               cc_manual_to_html -o ../../doc_html/basic/TriangulationDS_3  \
	         TriangulationDS_3/main.tex
	LATEX_CONV_INPUTS=.:Triangulation_3:TriangulationDS_3:Triangulation_3_ref:TDS2_ref\
              cc_manual_to_html  -o ../../doc_html/basic/TriangulationDS_3_ref \
	         TriangulationDS_3_ref/main.tex
	cp Triangulation_3/*.gif ../../doc_html/basic/Triangulation_3
	cp Triangulation_3_ref/*.gif ../../doc_html/basic/Triangulation_3_ref
	cp TriangulationDS_3/*.gif ../../doc_html/basic/TriangulationDS_3
	cp TriangulationDS_3_ref/*.gif ../../doc_html/basic/TriangulationDS_3_ref

#
#
#
#
