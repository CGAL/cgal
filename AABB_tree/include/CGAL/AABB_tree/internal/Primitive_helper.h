// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Sebastien Loriot

#ifndef CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER
#define CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER

#include <CGAL/license/AABB_tree.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian/Aff_transformation_3.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/AABB_tree/internal/Has_nested_type_Shared_data.h>
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


template<class Kernel>
Bbox_3 compute_transformed_bbox(const CGAL::Aff_transformation_3<Kernel>& at, const Bbox_3& bbox, bool has_rotation)
{
  typedef Simple_cartesian<Interval_nt<false>> AK;
  typedef Cartesian_converter<Kernel, AK>    C2F;
  C2F c2f;
  CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_UPWARD);

  AK::Aff_transformation_3 a_at = c2f(at);
  AK::FT xtrm[6] = { c2f((bbox.min)(0)), c2f((bbox.max)(0)),
                     c2f((bbox.min)(1)), c2f((bbox.max)(1)),
                     c2f((bbox.min)(2)), c2f((bbox.max)(2)) };

  if(!has_rotation){
    AK::Point_3 ps[2];
    ps[0] = a_at( AK::Point_3(xtrm[0], xtrm[2], xtrm[4]) );
    ps[1] = a_at( AK::Point_3(xtrm[1], xtrm[3], xtrm[5]) );

    return bbox_3(ps, ps+2);
  }

  AK::Point_3 ps[8];
  ps[0] = a_at( AK::Point_3(xtrm[0], xtrm[2], xtrm[4]) );
  ps[1] = a_at( AK::Point_3(xtrm[0], xtrm[2], xtrm[5]) );
  ps[2] = a_at( AK::Point_3(xtrm[0], xtrm[3], xtrm[4]) );
  ps[3] = a_at( AK::Point_3(xtrm[0], xtrm[3], xtrm[5]) );

  ps[4] = a_at( AK::Point_3(xtrm[1], xtrm[2], xtrm[4]) );
  ps[5] = a_at( AK::Point_3(xtrm[1], xtrm[2], xtrm[5]) );
  ps[6] = a_at( AK::Point_3(xtrm[1], xtrm[3], xtrm[4]) );
  ps[7] = a_at( AK::Point_3(xtrm[1], xtrm[3], xtrm[5]) );

  return bbox_3(ps, ps+8);
}

template<class Kernel>
Bbox_2 compute_transformed_bbox(const CGAL::Aff_transformation_2<Kernel>& at, const Bbox_2& bbox, bool has_rotation)
{
  typedef Simple_cartesian<Interval_nt<false>> AK;
  typedef Cartesian_converter<Kernel, AK>    C2F;
  C2F c2f;


  AK::Aff_transformation_2 a_at = c2f(at);
  AK::FT xtrm[4] = { c2f((bbox.min)(0)), c2f((bbox.max)(0)),
                     c2f((bbox.min)(1)), c2f((bbox.max)(1)) };

  if(!has_rotation){
    AK::Point_2 ps[2];
    ps[0] = a_at( AK::Point_2(xtrm[0], xtrm[2]) );
    ps[1] = a_at( AK::Point_2(xtrm[1], xtrm[3]) );

    return bbox_2(ps, ps+2);
  }

  AK::Point_2 ps[4];
  ps[0] = a_at( AK::Point_2(xtrm[0], xtrm[2]) );
  ps[1] = a_at( AK::Point_2(xtrm[0], xtrm[3]) );
  ps[2] = a_at( AK::Point_2(xtrm[1], xtrm[2]) );
  ps[3] = a_at( AK::Point_2(xtrm[1], xtrm[3]) );

  return bbox_2(ps, ps+4);
}

} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER
