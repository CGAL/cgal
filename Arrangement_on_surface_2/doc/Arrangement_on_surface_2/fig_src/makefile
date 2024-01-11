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
  vertical_decomposition.pdf \
  compare_x_on_boundary.pdf \
  compare_x_near_boundary.pdf \
  compare_y_near_boundary.pdf \
  identified_clusters.pdf \
  generic_polyline.pdf \
  well_oriented_polyline.pdf \
  circular_arcs.pdf \
  conic_multiplicities.pdf \
  rational_function_singular.pdf \
  consolidated_curve_data.pdf \
  generic_curve_data.pdf \
  observer.pdf \
  dcel_extension.pdf \
  overlay.pdf \
  overlay_unbounded.pdf \
  overlay_color.pdf \
  bgl_primal_adapter.pdf \
  bgl_dual_adapter.pdf \
  unb_asymptote.pdf \
  min_area_tri.pdf \
  jeep.pdf

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
