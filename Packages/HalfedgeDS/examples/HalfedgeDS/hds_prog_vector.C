// hds_prog_vector.C              
// -------------------------------------
#include <CGAL/HalfedgeDS_items.h>
#include <CGAL/HalfedgeDS_using_vector.h>
#include <CGAL/HalfedgeDS_decorator.h>

struct Traits { typedef int Point; };
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
 typedef CGAL::HalfedgeDS_using_vector     <Traits,CGAL::HalfedgeDS_items> HDS;
#else
 typedef CGAL::HalfedgeDS_using_vector::HDS<Traits,CGAL::HalfedgeDS_items> HDS;
#endif
typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;

int main() {
    HDS hds(1,2,2);
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
