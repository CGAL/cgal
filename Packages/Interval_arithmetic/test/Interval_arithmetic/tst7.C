#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Point_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Vector_2.h>

#include <vector>

#include <CGAL/Interval_arithmetic/_FPU.h>

CGAL_BEGIN_NAMESPACE
template <class FT>
inline
FT
square  (const VectorC2<FT> & v)
{ return v*v; }
CGAL_END_NAMESPACE

void pipo()
{
  CGAL::PointC2<double> p(0,1), q(1,1);
  double b = CGAL::square(p-q);
  std::cout << b << std::endl;
}

template < class NT>
void foo( CGAL::Ray_2< CGAL::Cartesian<NT> >& r) {
    typedef CGAL::Cartesian<NT>   R;
    typedef CGAL::Point_2<R>      Point;
    typedef CGAL::Direction_2<R>  Direction;
    typedef CGAL::Ray_2<R>        Ray;

    Point      pt;
    Direction  dir;
    r = Ray(pt,dir);
}

typedef CGAL::Cartesian<double> REP;
typedef CGAL::Ray_2<REP>        Ray;

int main ()
{
pipo();
    int i;
    for (i=0; i<200000000; i++)
        CGAL::FPU_set_cw(CGAL::FPU_cw_up);
    Ray var;
    foo( var);
    return 0;
}
