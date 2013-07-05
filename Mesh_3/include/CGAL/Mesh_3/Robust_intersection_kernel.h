// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_3/Robust_intersection_kernel.h $
// $Id: Robust_intersection_kernel.h 67573 2012-02-02 14:54:51Z lrineau $
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_ROBUST_INTERSECTION_KERNEL_3_H
#define CGAL_MESH_3_ROBUST_INTERSECTION_KERNEL_3_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>


namespace CGAL {

namespace Mesh_3 {

template < typename K_ >
class Robust_intersection_for_kernel_3
{
public:
  typedef typename K_::Line_3                         Line_3;
  typedef typename K_::Plane_3                        Plane_3;
  
  typedef Object                                      result_type;

  typedef Robust_intersection_for_kernel_3<K_>        Self;

  typedef Exact_predicates_exact_constructions_kernel EK;
  // typedef Simple_cartesian<Gmpq>  EK;
  typedef Cartesian_converter<typename K_::Kernel, EK>  To_exact;
  typedef Cartesian_converter<EK, typename K_::Kernel>  Back_from_exact;

  // forward to the base
  template<class T1, class T2>
  Object operator() (const T1& t, const T2& s) const
  {
    return CGAL::CommonKernelFunctors::Intersect_3<K_>()(t, s);
  }

  // exact computation only for Intersect_3()(Line_3, Plane_3)
  Object operator() (const Line_3& line, const Plane_3& plane) const
  {
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Intersect_3 exact_intersection = EK().intersect_3_object();
    Object object = exact_intersection(to_exact(line), to_exact(plane));
    // std::exit(1);
    if ( const EK::Point_3* p = object_cast<EK::Point_3>(&object) )
      return make_object(back_from_exact(*p));
    else if ( const EK::Segment_3* seg = object_cast<EK::Segment_3>(&object) )
      return make_object(back_from_exact(*seg));
    else 
      return Object();
  }
}; // end template Robust_intersection_for_kernel_3


template <typename K_base, typename Kernel>
struct Robust_intersection_kernel_base
  : public K_base::template Base<Kernel>::Type
{
  
  // template < typename Kernel2 >
  // struct Base {
  //   typedef typename K_base::template Base<Kernel2> K2;
  //   typedef Robust_intersection_kernel_base<K2>  Type;
  // };

  typedef Robust_intersection_for_kernel_3<Kernel> Intersect_3;

  Intersect_3
  intersect_3_object() const
  {
    return Intersect_3();
  }
}; // end template Robust_intersection_kernel_base

template <typename K>
struct Robust_intersection_kernel
  : Type_equality_wrapper<Robust_intersection_kernel_base<K,
                                                          Robust_intersection_kernel<K> >,
                          Robust_intersection_kernel<K> >
{};



} // end namespace Mesh_3
  
} //namespace CGAL

#endif // CGAL_MESH_3_ROBUST_INTERSECTION_KERNEL_3_H
