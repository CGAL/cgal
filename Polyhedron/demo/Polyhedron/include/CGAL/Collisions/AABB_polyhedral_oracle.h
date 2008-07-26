#ifndef CGAL_AABB_POLYHEDRAL_ORACLE_H
#define CGAL_AABB_POLYHEDRAL_ORACLE_H

#include <boost/static_warning.hpp>
#include <utility>
#include <CGAL/iterator.h>

namespace CGAL {

template <class Surface, class Kernel>
class AABB_polyhedral_oracle
{
public:
  typedef Surface Surface_3;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Ray_3 Ray;
  typedef typename Kernel::Line_3 Line;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename AABB_polyhedral_oracle<Surface,Kernel> Self;
  typedef typename Self Surface_mesher_traits_3;

public:

  // Surface constructor
  AABB_polyhedral_oracle()
  {
  }

  class Intersect_3;

  friend class Intersect_3;

  class Intersect_3 {
    const Self& self;
  public:
    Intersect_3(const Self& self) : self(self)
    {
    }

    Object operator()(const Surface_3& surface, const Segment_3& s) const
    {
      return Object();
    }
    
    Object operator()(const Surface_3& surface, const Ray_3& r) const
    {
      return Object();
    }
      
    Object operator()(const Surface_3& surface, const Line_3& l) const
    {
      return Object();
    }
  };

  Intersect_3 intersect_3_object() const
  {
    return Intersect_3(*this);
  }

  class Construct_initial_points;

  friend class Construct_initial_points;

  class Construct_initial_points
  {
    const Self& self;
  public:
    Construct_initial_points(const Self& self) : self(self)
    {
    }

    template <typename OutputIteratorPoints>
    OutputIteratorPoints operator() (const Surface_3& surface, 
                                     OutputIteratorPoints out, 
                                     int n) const 
    {
      // *out++= p;
      return out;
    }
  };

  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  template <class P>
  bool is_in_volume(const Surface_3& surface, const P& p)
  {
    return true;
  }
}; // end class AABB_polyhedral_oracle

} // end namespace CGAL

#endif // CGAL_AABB_POLYHEDRAL_ORACLE_H
