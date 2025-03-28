#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

#include <iostream>

typedef CGAL::Cartesian<CGAL::Exact_rational>                           Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<Kernel>        ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                 Traits;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                   Domain;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>           Factory;

typedef typename Traits::FT                                             FT;
typedef typename Traits::Hyperbolic_point_2                             Point;
typedef typename Traits::Complex                                        Complex;

int main()
{
  Factory factory;
  Domain domain = factory.make_hyperbolic_fundamental_domain_g2(3459);

  assert(domain.size()==8);
  for (std::size_t k=0; k<8; ++k) {
    assert(domain.paired_side(k)==(k+4)%8);
  }

  return 0;
}
