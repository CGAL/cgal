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

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Bbox_3.h>

#include <cmath>

// Test intersection of two bbox according to two transformation without rotation
template <typename Kernel>
bool do_intersect_transformed_BB(const CGAL::Bbox_3& b1,
                                 const CGAL::Bbox_3& b2,
                                 const CGAL::Aff_transformation_3<Kernel>& t1,
                                 const CGAL::Aff_transformation_3<Kernel>& t2)
{
  typedef Simple_cartesian<Interval_nt_advanced> AK;
  typedef Cartesian_converter<Kernel, AK>    C2F;
  C2F c2f;

  AK::Aff_transformation_3 a_t1 = c2f(t1);
  AK::FT xtrm1[6] = {c2f((b1.min)(0)), c2f((b1.max)(0)),
                     c2f((b1.min)(1)), c2f((b1.max)(1)),
                     c2f((b1.min)(2)), c2f((b1.max)(2)) };

  AK::Aff_transformation_3 a_t2 = c2f(t2);
  AK::FT xtrm2[6] = {c2f((b2.min)(0)), c2f((b2.max)(0)),
                     c2f((b2.min)(1)), c2f((b2.max)(1)),
                     c2f((b2.min)(2)), c2f((b2.max)(2)) };

  AK::Point_3 ps[4];
  ps[0] = a_t1( AK::Point_3(xtrm1[0], xtrm1[2], xtrm1[4]) );
  ps[1] = a_t1( AK::Point_3(xtrm1[1], xtrm1[3], xtrm1[5]) );
  ps[2] = a_t2( AK::Point_3(xtrm2[0], xtrm2[2], xtrm2[4]) );
  ps[3] = a_t2( AK::Point_3(xtrm2[1], xtrm2[3], xtrm2[5]) );

  return do_overlap(bbox_3(ps, ps+2), bbox_3(ps+2, ps+4));
}

// Tests if two oriented bounding boxes (OBBs) intersect using the Separating Axis Theorem (SAT).
// Tests separation along the 6 principal axes (3 from each OBB).
// Note: Does not test the 9 cross-product axes for efficiency; thus false positives may be returned.
// TODO: This could be accelerated using static filters instead of Intervals
template <typename Kernel>
bool do_intersect_OBB(const Bbox_3& b1,
                      const Bbox_3& b2,
                      const CGAL::Aff_transformation_3<Kernel>& t1,
                      const CGAL::Aff_transformation_3<Kernel>& t2)
{
  using AK = Simple_cartesian<Interval_nt_advanced>;
  using C2F = Cartesian_converter<Kernel, AK>;
  C2F c2f;

  using FT     = AK::FT;
  using Point  = AK::Point_3;
  using Vector = AK::Vector_3;

  AK::Aff_transformation_3 a_t1 = c2f(t1);
  AK::Aff_transformation_3 a_t2 = c2f(t2);

  // The center of each box.
  Point p1 = a_t1.transform( Point((FT(b1.xmin())+FT(b1.xmax()))/2, (FT(b1.ymin())+FT(b1.ymax()))/2, (FT(b1.zmin())+FT(b1.zmax()))/2) );
  Point p2 = a_t2.transform( Point((FT(b2.xmin())+FT(b2.xmax()))/2, (FT(b2.ymin())+FT(b2.ymax()))/2, (FT(b2.zmin())+FT(b2.zmax()))/2) );

  // Half width vectors of each box.
  const Vector A[3] = {
    a_t1.transform(Vector(b1.x_span()/2, 0, 0)),
    a_t1.transform(Vector(0, b1.y_span()/2, 0)),
    a_t1.transform(Vector(0, 0, b1.z_span()/2))
  };

  const Vector B[3] = {
    a_t2.transform(Vector(b2.x_span()/2, 0, 0)),
    a_t2.transform(Vector(0, b2.y_span()/2, 0)),
    a_t2.transform(Vector(0, 0, b2.z_span()/2))
  };

  const Vector dir = p2 - p1;

  // Test separation on the 3 axes of A.
  for(int i=0; i<3; ++i)
  {
    // Project the center and half-width along each axis and compare the distance between centers with the sum of half-width projections.
    const Vector& axis = A[i];
    const FT dist = CGAL::abs(dir * axis);

    const FT ra = axis.squared_length();
    const FT rb = CGAL::abs(B[0]*axis) + CGAL::abs(B[1]*axis) + CGAL::abs(B[2]*axis);
    if(dist > ra + rb) return false;
  }

  // Test separation on the 3 axes of B.
  for(int i = 0; i < 3; ++i)
  {
    const Vector& axis = B[i];
    const FT dist = CGAL::abs(dir * axis);

    const FT rb = axis.squared_length();
    const FT ra = CGAL::abs(A[0]*axis) + CGAL::abs(A[1]*axis) + CGAL::abs(A[2]*axis);
    if(dist > ra + rb) return false;
  }

  // No separating axis among the 6 face normals.
  return true;
}

template <typename Kernel>
bool do_intersect_transformed_BB(const CGAL::Bbox_2& b1,
                                 const CGAL::Bbox_2& b2,
                                 const CGAL::Aff_transformation_2<Kernel>& t1,
                                 const CGAL::Aff_transformation_2<Kernel>& t2)
{
  typedef Simple_cartesian<Interval_nt_advanced> AK;
  typedef Cartesian_converter<Kernel, AK>    C2F;
  C2F c2f;

  AK::Aff_transformation_2 a_t1 = c2f(t1);
  AK::FT xtrm1[4] = {c2f((b1.min)(0)), c2f((b1.max)(0)),
                     c2f((b1.min)(1)), c2f((b1.max)(1)) };

  AK::Aff_transformation_2 a_t2 = c2f(t2);
  AK::FT xtrm2[4] = {c2f((b2.min)(0)), c2f((b2.max)(0)),
                     c2f((b2.min)(1)), c2f((b2.max)(1)) };

  AK::Point_2 ps[4];
  ps[0] = a_t1( AK::Point_2(xtrm1[0], xtrm1[2]) );
  ps[1] = a_t1( AK::Point_2(xtrm1[1], xtrm1[3]) );
  ps[2] = a_t2( AK::Point_2(xtrm2[0], xtrm2[2]) );
  ps[3] = a_t2( AK::Point_2(xtrm2[1], xtrm2[3]) );

  return do_overlap(bbox_2(ps, ps+2), bbox_2(ps+2, ps+4));
}

// Tests if two oriented bounding boxes (OBBs) intersect using the Separating Axis Theorem (SAT).
// Tests separation along the 4 principal axes (2 from each OBB).
template <typename Kernel>
bool do_intersect_OBB(const Bbox_2& b1,
                      const Bbox_2& b2,
                      const CGAL::Aff_transformation_2<Kernel>& t1,
                      const CGAL::Aff_transformation_2<Kernel>& t2)
{
  using AK = Simple_cartesian<Interval_nt_advanced>;
  using C2F = Cartesian_converter<Kernel, AK>;
  C2F c2f;

  using FT     = AK::FT;
  using Point  = AK::Point_2;
  using Vector = AK::Vector_2;

  AK::Aff_transformation_2 a_t1 = c2f(t1);
  AK::Aff_transformation_2 a_t2 = c2f(t2);

  // The center of each box.
  Point p1 = a_t1.transform( Point((FT(b1.xmin())+FT(b1.xmax()))/2, (FT(b1.ymin())+FT(b1.ymax()))/2) );
  Point p2 = a_t2.transform( Point((FT(b2.xmin())+FT(b2.xmax()))/2, (FT(b2.ymin())+FT(b2.ymax()))/2) );

  // Half width vectors of each box.
  const Vector A[2] = {
    a_t1.transform(Vector(b1.x_span()/2, 0)),
    a_t1.transform(Vector(0, b1.y_span()/2))
  };

  const Vector B[2] = {
    a_t2.transform(Vector(b2.x_span()/2, 0)),
    a_t2.transform(Vector(0, b2.y_span()/2))
  };

  const Vector dir = p2 - p1;

  // Test separation on the 2 axes of A.
  for(int i=0; i<2; ++i)
  {
    // Project the center and half-width along each axis and compare the distance between centers with the sum of half-width projections.
    const Vector& axis = A[i];
    const FT dist = CGAL::abs(dir * axis);

    const FT ra = axis.squared_length();
    const FT rb = CGAL::abs(B[0]*axis) + CGAL::abs(B[1]*axis);
    if(dist > ra + rb) return false;
  }

  // Test separation on the 2 axes of B.
  for(int i = 0; i < 2; ++i)
  {
    const Vector& axis = B[i];
    const FT dist = CGAL::abs(dir * axis);

    const FT rb = axis.squared_length();
    const FT ra = CGAL::abs(A[0]*axis) + CGAL::abs(A[1]*axis);
    if(dist > ra + rb) return false;
  }

  // No separating axis among the 4 face normals.
  return true;
}

template<class Kernel>
Bbox_3 compute_transformed_bbox(const CGAL::Aff_transformation_3<Kernel>& at, const Bbox_3& bbox, bool has_rotation)
{
  typedef Simple_cartesian<Interval_nt_advanced> AK;
  typedef Cartesian_converter<Kernel, AK>    C2F;
  C2F c2f;

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
Bbox_3 compute_transformed_bbox(const CGAL::Aff_transformation_3<Kernel>& at, const Bbox_3& bbox)
{
  return compute_transformed_bbox(at, bbox, at.has_rotation());
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

template<class Kernel>
Bbox_2 compute_transformed_bbox(const CGAL::Aff_transformation_2<Kernel>& at, const Bbox_2& bbox)
{
  return compute_transformed_bbox(at, bbox, at.has_rotation());
}

} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_AABB_TREE_PRIMITIVE_HELPER
