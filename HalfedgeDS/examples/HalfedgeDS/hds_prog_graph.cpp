#include <CGAL/HalfedgeDS_min_items.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>

// no traits needed, argument can be arbitrary dummy.
typedef CGAL::HalfedgeDS_default<int, CGAL::HalfedgeDS_min_items> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
