// Copyright (c) 2011 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H
#define CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H

#include <CGAL/property_map.h>
#include <boost/type_traits/is_const.hpp>
#include <boost/mpl/if.hpp>

namespace CGAL{

//property map
template <class Polyhedron> 
struct Triangle_from_facet_property_map{
  //classical typedefs
  typedef typename boost::mpl::if_<
    typename boost::is_const<Polyhedron>::type,
    const typename Polyhedron::Facet,
    typename Polyhedron::Facet >::type key_type;
  typedef typename Polyhedron::Traits::Kernel::Triangle_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
};
//get function for property map
template <class Polyhedron>
inline
typename Polyhedron::Traits::Kernel::Triangle_3
get(const Triangle_from_facet_property_map<Polyhedron>&,
    typename Triangle_from_facet_property_map<Polyhedron>::key_type& f)
{
  typedef typename Polyhedron::Traits::Kernel Kernel;
  CGAL_precondition(f.halfedge() == f.halfedge()->next()->next()->next());
  const typename Kernel::Point_3& a = f.halfedge()->vertex()->point();
  const typename Kernel::Point_3& b = f.halfedge()->next()->vertex()->point();
  const typename Kernel::Point_3& c = f.halfedge()->next()->next()->vertex()->point();  
  return typename Kernel::Triangle_3(a,b,c);
}


template <class Polyhedron> 
struct Segment_from_halfedge_property_map{
  //classical typedefs
  typedef typename boost::mpl::if_<
    typename boost::is_const<Polyhedron>::type,
    const typename Polyhedron::Halfedge,
    typename Polyhedron::Halfedge >::type key_type;
  typedef typename Polyhedron::Traits::Kernel::Segment_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
};
//get function for property map
template <class Polyhedron>
inline
typename Polyhedron::Traits::Kernel::Segment_3
get(const Segment_from_halfedge_property_map<Polyhedron>&,
    typename Segment_from_halfedge_property_map<Polyhedron>::key_type& h)
{
  typedef typename Polyhedron::Traits::Kernel Kernel;
  const typename Kernel::Point_3& a = h.vertex()->point();
  const typename Kernel::Point_3& b = h.opposite()->vertex()->point();
  return typename Kernel::Segment_3(a,b);
}

} //namespace CGAL

#endif //CGAL_POLYHEDRON_SIMPLEX_PROPERTY_MAP_H
