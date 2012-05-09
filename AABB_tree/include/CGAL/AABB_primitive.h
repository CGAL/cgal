// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/branches/features/AABB_tree-one_primitive_per_object-sloriot/AABB_tree/include/CGAL/AABB_triangle_primitive.h $
// $Id: AABB_triangle_primitive.h 68970 2012-05-04 14:59:14Z sloriot $
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_TRIANGLE_PRIMITIVE_H_
#define CGAL_AABB_TRIANGLE_PRIMITIVE_H_

#include <CGAL/internal/AABB_tree/Primitive_caching.h>
#include <CGAL/property_map.h>
#include <CGAL/Default.h>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/mpl/if.hpp>

namespace CGAL {

template < class Iterator,
           class ObjectPropertyMap,
           class PointPropertyMap,
           bool cache_primitive=false >
class AABB_primitive :
  public internal::Primitive_caching< Iterator, TrianglePropertyMap, cache_primitive >
{
  // types
  typedef internal::Primitive_caching<Iterator,TrianglePropertyMap,cache_primitive> Primitive_base;
public:
  
  typedef typename boost::property_traits< ObjectPropertyMap >::value_type Datum; //datum type
  typedef typename boost::property_traits< PointPropertyMap  >::value_type Point; //point type
  typedef Iterator Id; // Id type

private:
  Id m_it;
  PointPropertyMap m_ppmap
public:
  // constructors
  AABB_primitive() {}
  AABB_primitive(Id it,ObjectPropertyMap t_pmap=ObjectPropertyMap(), PointPropertyMap p_pmap=PointPropertyMap())
          : m_ppmap(p_pmap), m_it(it)
  {
    this->set_primitive(it,t_pmap);
  }
public:
  Id& id() { return m_it; }
  const Id& id() const { return m_it; }
  
  typename Primitive_base::result_type datum() const {
    return this->get_primitive(m_it);
  }

  /// Returns a point on the primitive
  typename boost::property_traits< PointPropertyMap  >::reference 
  reference_point() const {
    return get(m_ppmap,*m_it);
  }
};

}  // end namespace CGAL


#endif // CGAL_AABB_TRIANGLE_PRIMITIVE_H_

