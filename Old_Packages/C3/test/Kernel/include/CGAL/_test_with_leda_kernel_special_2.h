#line 390 "test_kernel_programs.fw"
#include <CGAL/Point_2.h>
#include <LEDA/rat_point.h>

class CGAL_use_rat_leda_kernel
{
  public:
    typedef  leda_integer      RT;
    typedef  leda_rational     FT;

    typedef  leda_rat_point    Point_2;
    typedef  leda_rat_vector   Vector_2;
};


CGAL_TEMPLATE_NULL
class CGAL_Point_2<CGAL_use_rat_leda_kernel> : public leda_rat_point 
{
  public:

    typedef CGAL_use_rat_leda_kernel   R;
    typedef typename R::RT             RT;
    typedef typename R::FT             FT;

    CGAL_Point_2<CGAL_use_rat_leda_kernel>()
    {}
  
    CGAL_Point_2<CGAL_use_rat_leda_kernel>(const CGAL_Origin &o)
    {}

    CGAL_Point_2<CGAL_use_rat_leda_kernel>(const leda_rat_point& p)
      : leda_rat_point(p)
    {}

    CGAL_Point_2<CGAL_use_rat_leda_kernel>(const RT &hx, const RT &hy)
      : leda_rat_point(hx, hy)
    {}

    CGAL_Point_2<CGAL_use_rat_leda_kernel>(const RT &hx, 
                                           const RT &hy, 
                                           const RT &hw)
      : leda_rat_point(hx, hy, hw)
    {}
};
