#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/HalfedgeDS_decorator.h>

struct Traits { typedef int Point_2; };
typedef CGAL::HalfedgeDS_vector< Traits, CGAL::HalfedgeDS_items_2> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;

int main() {
    HDS hds(1,2,2);
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
