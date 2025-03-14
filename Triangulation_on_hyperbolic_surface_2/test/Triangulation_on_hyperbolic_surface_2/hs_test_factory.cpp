#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

#include <iostream>
#include <vector>

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

  std::vector<Point> vertices;
  Point z0 = Point(FT("2057/2500"),FT("0"));
  Point z1 = Point(FT("6183/10000"),FT("2369/5000"));
  Point z2 = Point(FT("93/625"),FT("9299/10000"));
  Point z3 = Point(FT("-129263137397146717229370666475391966694612741374/530219238243800202784978257516468117328515435625"),FT("443541103604461956104066082668282544568945930056/530219238243800202784978257516468117328515435625"));
  Point z4 = Point(FT("-2057/2500"),FT("0"));
  Point z5 = Point(FT("-6183/10000"),FT("-2369/5000"));
  Point z6 = Point(FT("-93/625"),FT("-9299/10000"));
  Point z7 = Point(FT("129263137397146717229370666475391966694612741374/530219238243800202784978257516468117328515435625"),FT("-443541103604461956104066082668282544568945930056/530219238243800202784978257516468117328515435625"));
  vertices.push_back(z0);
  vertices.push_back(z1);
  vertices.push_back(z2);
  vertices.push_back(z3);
  vertices.push_back(z4);
  vertices.push_back(z5);
  vertices.push_back(z6);
  vertices.push_back(z7);

  assert(domain.size()==8);
  for (std::size_t k=0; k<8; ++k) {
    assert(domain.vertex(k)==vertices[k]);
    assert(domain.paired_side(k)==(k+4)%8);
  }

  return 0;
}
