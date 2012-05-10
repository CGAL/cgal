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
// $URL$
// $Id$
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_PRIMITIVE_H
#define CGAL_AABB_PRIMITIVE_H

#include <CGAL/internal/AABB_tree/Primitive_caching.h>

namespace CGAL {

template < class Id_,
           class ObjectPropertyMap,
           class PointPropertyMap,
           class cache_datum=Tag_false >
class AABB_primitive :
  public internal::Primitive_caching< Id_, ObjectPropertyMap, cache_datum >
{
  // types
  typedef internal::Primitive_caching<Id_, ObjectPropertyMap, cache_datum> Primitive_base;
public:
  
  typedef typename boost::property_traits< ObjectPropertyMap >::value_type Datum; //datum type
  typedef typename boost::property_traits< PointPropertyMap  >::value_type Point; //point type
  typedef Id_ Id; // Id type

private:
  PointPropertyMap m_ppmap;
  Id m_it;

public:
  // constructors
  AABB_primitive(Id it,ObjectPropertyMap t_pmap=ObjectPropertyMap(), PointPropertyMap p_pmap=PointPropertyMap())
          : Primitive_base(it,t_pmap),m_ppmap(p_pmap), m_it(it)
  {}
public:
  Id& id() { return m_it; }
  const Id& id() const { return m_it; }
  
  typename Primitive_base::result_type datum() const {
    return this->get_primitive(m_it);
  }

  /// Returns a point on the primitive
  typename boost::property_traits< PointPropertyMap  >::reference 
  reference_point() const {
    return get(m_ppmap, m_it);
  }
};

}  // end namespace CGAL


#endif // CGAL_AABB_PRIMITIVE_H

