all : doc_ps doc_html

LOCALINPUTS=.:../../examples:Alpha_shapes_3:Alpha_shapes_3_ref:

doc_ps :
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} latex wrapper
	bibtex wrapper
	makeindex wrapper
	index_fix wrapper.ind
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} latex wrapper
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} dvips wrapper -o

doc_html :
	rm -rf  ../../doc_html/basic/
	mkdir -p ../../doc_html/basic/Alpha_shapes_3
	mkdir -p ../../doc_html/basic/Alpha_shapes_3_ref
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
             cc_manual_to_html -o ../../doc_html/basic/Alpha_shapes_3\
	         Alpha_shapes_3/main.tex
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
              cc_manual_to_html  -o ../../doc_html/basic/Alpha_shapes_3_ref \
                 Alpha_shapes_3_ref/main.tex 
	cp Alpha_shapes_3/*.gif ../../doc_html/basic/Alpha_shapes_3


#
#
#
#
