// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Sebastien Loriot

#ifndef CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER
#define CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER

#include <CGAL/license/AABB_tree.h>


#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL{
namespace internal{

//for backward compatibility: if Datum_reference and Point_reference are not defined in the primitive
//(using auto would solve the pb)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Datum_reference,Datum_reference,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Point_reference,Point_reference,false)

template<class Primitive,bool has_nested_type=Has_nested_type_Datum_reference<Primitive>::value>
struct Datum_result_type{ typedef typename Primitive::Datum_reference type; };
template<class Primitive>
struct Datum_result_type<Primitive,false>{ typedef typename Primitive::Datum type; };
template<class Primitive,bool has_nested_type=Has_nested_type_Point_reference<Primitive>::value>
struct Point_result_type{ typedef typename Primitive::Point_reference type; };
template<class Primitive>
struct Point_result_type<Primitive,false>{ typedef typename Primitive::Point type; };


//helper controlling whether extra data should be stored in the AABB_tree traits class  
template <class AABBTraits, bool has_shared_data=Has_nested_type_Shared_data<typename AABBTraits::Primitive>::value>
struct Primitive_helper;

template <class AABBTraits>
struct Primitive_helper<AABBTraits,true>{
  typedef typename Datum_result_type<typename AABBTraits::Primitive>::type Datum_type;
  static Datum_type get_datum(const typename AABBTraits::Primitive& p,const AABBTraits& traits)
  {
    return p.datum(traits.shared_data());
  }
  typedef typename Point_result_type<typename AABBTraits::Primitive>::type Reference_point_type;
  static Reference_point_type get_reference_point(const typename AABBTraits::Primitive& p,const AABBTraits& traits) {
    return p.reference_point(traits.shared_data());
  }
};
  
template <class AABBTraits>
struct Primitive_helper<AABBTraits,false>{
  typedef typename Datum_result_type<typename AABBTraits::Primitive>::type Datum_type;
  static Datum_type get_datum(const typename AABBTraits::Primitive& p,const AABBTraits&) {return p.datum();}
  typedef typename Point_result_type<typename AABBTraits::Primitive>::type Reference_point_type;
  static Reference_point_type get_reference_point(const typename AABBTraits::Primitive& p,const AABBTraits&) {return p.reference_point();}
};

} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER
