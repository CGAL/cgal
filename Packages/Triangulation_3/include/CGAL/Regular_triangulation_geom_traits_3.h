#ifndef CGAL_REGULAR_TRIANGULATION_GEOM_TRAITS_3_H
#define CGAL_REGULAR_TRIANGULATION_GEOM_TRAITS_3_H

#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Weighted_point_3.h>


template < class R , class W >
class CGAL_Regular_triangulation_geom_traits_3
  : public CGAL_Triangulation_geom_traits_3<R>
{
public:
  typedef R                                     Rep;
  typedef W                                     Weight;
  typedef CGAL_Triangulation_geom_traits_3 <R>  Geom;
  typedef Geom::Point                           Bare_point;
  typedef Geom::Segment                         Segment;
  typedef Geom::Triangle                        Triangle;
  typedef Geom::Tetrahedron                     Tetrahedron;
  typedef CGAL_Weighted_point_3<Bare_point, W>  Weighted_point;
  typedef Weighted_point                        Point;

  CGAL_Regular_triangulation_geom_traits_3()
    : Geom()
    {}

  
  CGAL_Oriented_side power_test_3(const Weighted_point &c1,
				  const Weighted_point &c2,
				  const Weighted_point &c3,
				  const Weighted_point &c4,
				  const Weighted_point &test) const
    {
      return CGAL_power_test_3(c1,c2,c3,c4,test);
    }

  CGAL_Oriented_side power_test_3(const Weighted_point &c1,
				  const Weighted_point &c2,
				  const Weighted_point &c3,
				  const Weighted_point &test) const
    {
      return CGAL_power_test_3(c1,c2,c3,test);
    }

  CGAL_Oriented_side power_test_3(const Weighted_point &c1,
				  const Weighted_point &c2,
				  const Weighted_point &test) const
    {
      return CGAL_power_test_3(c1,c2,test);
    }
};

#endif

