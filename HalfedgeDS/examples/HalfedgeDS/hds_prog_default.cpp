#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <cassert>

struct Traits { typedef int Point_2; };
typedef CGAL::HalfedgeDS_default<Traits> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    assert( decorator.is_valid());
    return 0;
}
