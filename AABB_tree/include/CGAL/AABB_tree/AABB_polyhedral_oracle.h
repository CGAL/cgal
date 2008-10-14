// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Camille Wormser, Pierre Alliez

#ifndef CGAL_AABB_POLYHEDRAL_ORACLE_H
#define CGAL_AABB_POLYHEDRAL_ORACLE_H

#include <boost/static_warning.hpp>
#include <utility>
#include <CGAL/iterator.h>

#include "AABB_tree.h"

namespace CGAL {

template <class Polyhedron, class Kernel, class AABBTree_kernel>
class AABB_polyhedral_oracle : public Polyhedron
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;

  typedef AABB_polyhedral_oracle<Polyhedron,Kernel,AABBTree_kernel> Self;
  typedef Self Surface_mesher_traits_3;
  typedef Point_3 Intersection_point;
  typedef Self Surface_3;

  // AABB tree
  typedef AABB_tree<AABBTree_kernel,typename Polyhedron::Facet_handle,Polyhedron> Tree;
  typedef typename Tree::Point_with_input Point_with_facet_handle;
  typedef CGAL::Cartesian_converter<Kernel,AABBTree_kernel> Converter;
  typedef CGAL::Cartesian_converter<AABBTree_kernel,Kernel> BConverter;
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
      Converter convert;
      BConverter bconvert;
      Point_with_facet_handle pwh;
      if(surface.tree()->first_intersection(convert(segment),pwh))
	return make_object(bconvert(pwh.first));
      else
	return Object();
    }
    
    Object operator()(const Surface_3& surface, const Ray_3& ray) const
    {
      Converter convert;
      BConverter bconvert;
      Point_with_facet_handle pwh;
      if(surface.tree()->first_intersection(convert(ray),pwh))
	return make_object(bconvert(pwh.first));
      else
	return Object();
    }
      
    Object operator()(const Surface_3& surface, const Line_3& line) const
    {
      Converter convert;
      BConverter bconvert;
      Point_with_facet_handle pwh;
      if(surface.tree()->first_intersection(convert(line),pwh))
	return make_object(bconvert(pwh.first));
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
			std::cout << "AABB_polyhedral_oracle: construct initial point set not implemented" << std::endl;
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
		std::cout << "call is in volume: empty function" << std::endl;
    return true;
  }
}; // end class AABB_polyhedral_oracle

} // end namespace CGAL

#endif // CGAL_AABB_POLYHEDRAL_ORACLE_H
