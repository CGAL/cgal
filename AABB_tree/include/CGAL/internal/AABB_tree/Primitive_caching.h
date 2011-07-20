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
//******************************************************************************
// File Description :
//
//******************************************************************************


#ifndef CGAL_INTERNAL_AABB_TREE_PRIMITIVE_CACHING_H
#define CGAL_INTERNAL_AABB_TREE_PRIMITIVE_CACHING_H

namespace CGAL {
namespace internal{

  template <class Primitive,class Id,class PropertyMap,bool do_cache>
  struct Primitive_caching;
  
  template <class Primitive,class Id,class PropertyMap>
  struct Primitive_caching<Primitive,Id,PropertyMap,true>
  {
    typedef const Primitive& result_type;
    Primitive datum;
    
    void set_primitive(Id id,PropertyMap pmap){datum=get(pmap,*id);}
    result_type get_primitive(Id) const{
      return datum;
    }
  };

  template <class Primitive,class Id,class PropertyMap>
  struct Primitive_caching<Primitive,Id,PropertyMap,false>
  {
    typedef Primitive result_type;
    PropertyMap pmap_;
    
    void set_primitive(Id,PropertyMap pmap){pmap_=pmap;}
    result_type get_primitive(Id id) const{
      return get(pmap_,*id);
    }
  };


} } //namespace CGAL::internal

#endif
