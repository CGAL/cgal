# latex ws.tex
# dvips -E ws.dvi -o
# epstool --bbox --copy ws.ps ws.eps
# ps2pdf -dEPSCrop ws.eps ws.pdf

TARGETS = \
  simple_arr.pdf \
  edge_manipulation.pdf \
  special_edge_insertion.pdf \
  isolated_vertices.pdf \
  edge_insertion.pdf \
  triangle.pdf \
  pred_around_vertex.pdf \
  arr_segs.pdf \
  insert.pdf \
  point_location.pdf \
  incremental_insertion.pdf \
  global_insertion.pdf \
  global_removal.pdf \
  unbounded_non_intersecting.pdf \
  unb_dcel.pdf \
  polylines.pdf \
  circles.pdf \
  conics.pdf \
  rational_functions.pdf \
  unbounded_rational_functions.pdf \
  bezier_curves.pdf \
  algebraic_curves.pdf \
  algebraic_curves_decomposition.pdf \
  algebraic_segments.pdf \
  curve_history.pdf \
  edge_manipulation_curve_history.pdf \
  vertical_decomposition.pdf

.DEFAULT_GOAL = all
.PHONY: all clean

all: $(TARGETS)

%.dvi: %.tex
	latex $<

%.ps: %.dvi
	dvips -E $< -o

%.eps: %.ps
	epstool --bbox --copy $< $@

%.pdf: %.eps
	ps2pdf -dEPSCrop $< $@

clean:
	rm -f $(DOCUMENTS)
