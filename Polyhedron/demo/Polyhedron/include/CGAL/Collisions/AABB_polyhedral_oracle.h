#ifndef CGAL_AABB_POLYHEDRAL_ORACLE_H
#define CGAL_AABB_POLYHEDRAL_ORACLE_H

#include <boost/static_warning.hpp>
#include <utility>
#include <CGAL/iterator.h>

#include "AABB_tree.h"

namespace CGAL {

template <class Polyhedron, class Kernel>
class AABB_polyhedral_oracle : public Polyhedron
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;

  typedef typename AABB_polyhedral_oracle<Polyhedron,Kernel> Self;
  typedef typename Self Surface_mesher_traits_3;
  typedef typename Point_3 Intersection_point;
  typedef Self Surface_3;
 // Visitor visitor;


  // AABB tree
  typedef AABB_tree<Kernel,typename Polyhedron::Facet_handle,Polyhedron> Tree;
  typedef typename Tree::Point_with_input Point_with_facet_handle;
  Tree *m_pTree;

public:
  Tree* tree() const { return m_pTree; }

public:
  // Surface constructor
  AABB_polyhedral_oracle()
  {
    m_pTree = NULL;
  }
  AABB_polyhedral_oracle(Tree *pTree)
  {
    m_pTree = pTree;
  }
  AABB_polyhedral_oracle(const AABB_polyhedral_oracle& oracle)
  {
    m_pTree = oracle.tree();
  }


  class Intersect_3;

  friend class Intersect_3;

  class Intersect_3 {
    const Self& self;
  public:
    Intersect_3(const Self& self) : self(self)
    {
    }

    Object operator()(const Surface_3& surface, const Segment_3& segment) const
    {
      Point_with_facet_handle pwh;
      if(surface.tree()->furthest_intersection(segment,CGAL::ORIGIN,pwh))
	//                                       ^^^^^^^^^^^^ to fix
	return make_object(pwh.first);
      else
	return Object();
    }
    
    Object operator()(const Surface_3& surface, const Ray_3& ray) const
    {
      Point_with_facet_handle pwh;
      if(surface.tree()->furthest_intersection(ray,CGAL::ORIGIN,pwh))
	//                                       ^^^^^^^^^^^^ to fix
	return make_object(pwh.first);
      else
	return Object();
    }
      
    Object operator()(const Surface_3& surface, const Line_3& line) const
    {
      Point_with_facet_handle pwh;
      if(surface.tree()->furthest_intersection(line,CGAL::ORIGIN,pwh))
	//                                       ^^^^^^^^^^^^ to fix
	return make_object(pwh.first);
      else
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
      // TODO (with visitor)
      // std::cout << "construct initial point set" << std::endl;
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
    std::cout << "DEBUG: call is in volume" << std::endl;
    return true;
  }
}; // end class AABB_polyhedral_oracle

} // end namespace CGAL

#endif // CGAL_AABB_POLYHEDRAL_ORACLE_H
