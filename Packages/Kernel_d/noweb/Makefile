INTERSECTIONOBJS = intersection_objectsHd.h intersection_objectsCd.h \
	intersection_objects_d.h 
FUNCTIONOBJS = function_objectsHd.h function_objectsCd.h 
KERNELTESTS = interface-test.C intersection-test.C afftrafo-test.C

KERNELSUB = ../include/CGAL/Kernel_d/
KERNELTOP = ../include/CGAL/
KERNELTESTDIR = ../test/Kernel_d/


tests : $(KERNELTESTS)
	cp $(KERNELTESTS) $(KERNELTESTDIR)

intobjs : $(INTERSECTIONOBJS)
	cp $(INTERSECTIONOBJS) $(KERNELSUB)

funcobjs : $(FUNCTIONOBJS)
	cp $(FUNCTIONOBJS) $(KERNELSUB)

$(KERNELTESTS) : Kernel_d-tests.lw 
	notangle -R$@ $< > $@

$(INTERSECTIONOBJS) : intersection_objects_d.lw
	notangle -R$@ $< > $@

function_objectsHd.h : function_objectsHd.lw
	notangle -R$@ $< > $@

function_objectsCd.h : function_objectsCd.lw
	notangle -R$@ $< > $@





