HYPERPLANE = HyperplaneHd.h HyperplaneHd.C
AFFTRAFO = Aff_transformationHd.h Aff_transformationHd.C
PREDS = 
INTERSECTIONOBJS = intersection_objectsHd.h
FUNCTIONOBJS = function_objectsHd.h
KERNELTESTS = interface-test.C intersection-test.C afftrafo-test.C

KERNELSUB = ../include/CGAL/Kernel_d/
KERNELTOP = ../include/CGAL/
KERNELTESTDIR = ../test/Kernel_d/


tests : $(KERNELTESTS)
	cp $(KERNELTESTS) $(KERNELTESTDIR)

hyperplane : $(HYPERPLANE)
	cp $(HYPERPLANE) $(KERNELSUB)

interobjs : $(INTERSECTIONOBJS)
	cp $(INTERSECTIONOBJS) $(KERNELSUB)

funcobjs : $(FUNCTIONOBJS)
	cp $(FUNCTIONOBJS) $(KERNELSUB)

afftrafo : $(AFFTRAFO)
	cp $(AFFTRAFO) $(KERNELSUB)

$(KERNELTESTS) : Kernel_d-tests.lw 
	notangle -R$@ $< > $@

$(INTERSECTIONOBJS) : IntersectionsHd.lw
	notangle -R$@ $< > $@

$(FUNCTIONOBJS) : function_objectsHD.lw
	notangle -R$@ $< > $@

$(HYPERPLANE) : HyperplaneHd.lw
	notangle -R$@ $< > $@

$(AFFTRAFO) : Aff_transformationHd.lw
	notangle -R$@ $< > $@





