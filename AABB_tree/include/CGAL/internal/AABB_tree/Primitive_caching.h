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

#include <CGAL/property_map.h>
#include <CGAL/tags.h>

#ifndef CGAL_INTERNAL_AABB_TREE_PRIMITIVE_CACHING_H
#define CGAL_INTERNAL_AABB_TREE_PRIMITIVE_CACHING_H

namespace CGAL {
namespace internal{

  template <class Id,class PropertyMap,class do_cache>
  struct Primitive_caching;
  
  template <class Id,class PropertyMap>
  struct Primitive_caching<Id,PropertyMap,Tag_true>
  {

    typedef typename boost::property_traits< PropertyMap >::value_type Primitive;
    typedef const Primitive& result_type;
    Primitive datum;
    
    Primitive_caching(Id id,PropertyMap pmap) {datum=get(pmap,id);}
    result_type get_primitive(Id) const{
      return datum;
    }
  };

  template <class Id,class PropertyMap>
  struct Primitive_caching<Id,PropertyMap,Tag_false>
  {
    typedef typename boost::property_traits< PropertyMap >::reference result_type;
    PropertyMap pmap_;
    
    Primitive_caching(Id,PropertyMap pmap){pmap_=pmap;}
    result_type get_primitive(Id id) const{
      return get(pmap_,id);
    }
  };


} } //namespace CGAL::internal

#endif
