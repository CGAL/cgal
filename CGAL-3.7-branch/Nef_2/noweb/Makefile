#---------------------------------------------------------------------#
# execution enforced by gmake -k ...
#---------------------------------------------------------------------#

MAININCL = ../include/CGAL
NEFINCL = $(MAININCL)/Nef_2
TESTDIR = ../test/Nef_2
DEMODIR = ../demo/Nef_2

# polynomials:
Polynomial-test.C Polynomial-demo.C : Polynomial.h
Polynomial.h : Polynomial.lw
	notangle -R$@ Polynomial.lw | specialize > $@

Polynomial-test.C Polynomial-demo.C : Polynomial.lw
	notangle -R$@ Polynomial.lw > $@

poly : Polynomial.h Polynomial-test.C 
	mv Polynomial.h $(NEFINCL)
	mv Polynomial-test.C $(TESTDIR)

Constrained_triang_traits.h Constrained_triang_anim.h \
	Constrained_triang-test.C : Constrained_triang.lw
	notangle -R$@ $< > $@

CTS = Constrained_triang_traits.h Constrained_triang_anim.h
ct : $(CTS)
	mv $(CTS) $(NEFINCL)

Extended_cartesian.h Extended_homogeneous.h Line_to_epoint.h \
	EPoint-demo.C EPoint-test.C : Simple_extended_kernel.lw

Filtered_extended_homogeneous.h : Filtered_extended_kernel.lw
	notangle -R$@ $< > $@

EKERNELS = Filtered_extended_homogeneous.h \
	Extended_cartesian.h Extended_homogeneous.h
ekernels : $(EKERNELS) EPoint-test.C
	mv $(EKERNELS) $(MAININCL)
	mv EPoint-test.C $(TESTDIR)

HDS_items.h : PM_decorator.lw
	notangle -R$@ $< > $@

items : HDS_items.h
	mv HDS_items.h $(NEFINCL)

Segment_overlay_traits.h \
leda_segment_overlay-test.c random_overlay-test.c \
leda_overlay_traits.h leda_overlay_visualization.h leda_segment_overlay-rt.c \
Affine_geometry.h HDS_decorator.h cgal_segment_overlay-test.C \
	: Segment_overlay.lw
	notangle -R$@ $< > $@

segint : Segment_overlay_traits.h
	mv Segment_overlay_traits.h $(NEFINCL)

PM_overlayer.h PM_overlayer-test.C PM_overlayer-demo.C : PM_overlayer.lw
	notangle -R$@ $< > $@

overlayer : PM_overlayer.h
	mv PM_overlayer.h $(NEFINCL)

PM_point_locator.h PM_persistent_PL.h PM_point_locator-test.C : \
	PM_point_locator.lw
	notangle -R$@ $< > $@

PLS = PM_point_locator.h PM_persistent_PL.h
pl : $(PLS)
	mv $(PLS) $(NEFINCL)


Nef_polyhedron_2.h PM_explorer.h Nef_polyhedron_2_Window_stream.h \
Nef_polyhedron_2-ps.C Nef_polyhedron_2-demo.C \
Nef_polyhedron_2-test.C Nef_polyhedron_2-rt.C : Nef_polyhedron_2.lw
	notangle -R$@ $< > $@

nef : Nef_polyhedron_2.h Nef_polyhedron_2-demo.C Nef_polyhedron_2-test.C
	mv Nef_polyhedron_2.h $(MAININCL)
	mv Nef_polyhedron_2-demo.C $(DEMODIR)
	mv Nef_polyhedron_2-test.C $(TESTDIR)

