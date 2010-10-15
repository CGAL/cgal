#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/HalfedgeDS_items_with_halfedge_cycle_2.h>

struct Traits { typedef int Point_2; };
typedef CGAL::HalfedgeDS_default<Traits,CGAL::HalfedgeDS_items_with_halfedge_cycle_2> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
