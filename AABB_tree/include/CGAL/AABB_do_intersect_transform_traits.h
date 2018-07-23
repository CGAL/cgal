
// Copyright (c) 2018 GeometryFactory (France).
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
//
// Author(s) : Maxime Gimeno
//

#ifndef CGAL_AABB_DO_INTERSECT_TRANSFORM_TRAITS_H
#define CGAL_AABB_DO_INTERSECT_TRANSFORM_TRAITS_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>
#include <CGAL/intersections.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/Filtered_predicate.h>

#include <CGAL/Aff_transformation_3.h>
#include <boost/optional.hpp>
#include <boost/bind.hpp>

/// \file AABB_do_intersect_transform_traits.h

//! \todo add protector
namespace CGAL {
// forward declaration
template< typename AABBTraits>
class AABB_tree;

namespace internal_AABB
{
template<class Kernel>
struct Transformed_datum_do_intersect_impl
{
  typedef bool result_type;
  template<class Query, class Datum_t>
  bool operator()(const Query& q,
                  const CGAL::Aff_transformation_3<Kernel>& transfo,
                  const Datum_t& pr) const
  {
    return Kernel().do_intersect_3_object()(q,
                                            pr.transform(transfo));
  }
};

template<class Kernel,
         class Has_filtered_predicates = typename Kernel::Has_filtered_predicates_tag /*Tag_false*/>
struct Transformed_datum_do_intersect
{
  template<class Query, class Datum_t>
  bool operator()(const Query& q,
                  const CGAL::Aff_transformation_3<Kernel>& transfo,
                  const Datum_t& pr) const
  {
    return Transformed_datum_do_intersect_impl<Kernel>()(q,transfo,pr);
  }
};

template<class Kernel>
struct Transformed_datum_do_intersect<Kernel, Tag_true>
{
  template<class Query, class Datum_t>
  bool operator()(const Query& q,
                  const CGAL::Aff_transformation_3<Kernel>& transfo,
                  const Datum_t& pr) const
  {
    typedef typename Kernel::Approximate_kernel     FK;
    typedef typename Kernel::Exact_kernel           EK;
    typedef typename Kernel::C2F                    C2F;
    typedef typename Kernel::C2E                    C2E;

    // filtered predicate
    CGAL::Filtered_predicate<
        internal_AABB::Transformed_datum_do_intersect_impl<EK>,
         internal_AABB::Transformed_datum_do_intersect_impl<FK>,
        C2E,
        C2F > filtered_do_intersect;

    return filtered_do_intersect(q, transfo, pr);
  }
};

}//end internal
template<typename BaseTraits,
         typename Kernel>
class AABB_do_intersect_transform_traits:
    public BaseTraits
{
  mutable Aff_transformation_3<Kernel> m_transfo;
public:

  //Constructor
  AABB_do_intersect_transform_traits(const Aff_transformation_3<Kernel>& transf = Aff_transformation_3<Kernel>(IDENTITY))
    :m_transfo(transf)
  {}
  // AABBTraits concept types
  typedef typename BaseTraits::Point_3 Point_3;
  typedef typename BaseTraits::Primitive Primitive;
  typedef typename BaseTraits::Bounding_box Bounding_box;
  //Intersections
  class Do_intersect
  {
    // TODO: possible optimization using Protector
    typedef Simple_cartesian<Interval_nt<> >             Approximate_kernel;
    typedef Cartesian_converter<Kernel, Approximate_kernel>    C2F;


    const AABB_do_intersect_transform_traits<BaseTraits, Kernel>& m_traits;
    C2F m_c2f;

    Bounding_box
    compute_transformed_bbox(const Bounding_box& bbox) const
    {
      Approximate_kernel::Aff_transformation_3 af = m_c2f( m_traits.transformation() );

      //TODO reuse the conversions
      typename Approximate_kernel::Point_3 ps[8];
      ps[0] = af( m_c2f( Point_3(bbox.min(0), bbox.min(1), bbox.min(2)) ) );
      ps[1] = af( m_c2f( Point_3(bbox.min(0), bbox.min(1), bbox.max(2)) ) );
      ps[2] = af( m_c2f( Point_3(bbox.min(0), bbox.max(1), bbox.min(2)) ) );
      ps[3] = af( m_c2f( Point_3(bbox.min(0), bbox.max(1), bbox.max(2)) ) );

      ps[4] = af( m_c2f( Point_3(bbox.max(0), bbox.min(1), bbox.min(2)) ) );
      ps[5] = af( m_c2f( Point_3(bbox.max(0), bbox.min(1), bbox.max(2)) ) );
      ps[6] = af( m_c2f( Point_3(bbox.max(0), bbox.max(1), bbox.min(2)) ) );
      ps[7] = af( m_c2f( Point_3(bbox.max(0), bbox.max(1), bbox.max(2)) ) );

      return bbox_3(ps, ps+8);
    }

  public:
    Do_intersect(const AABB_do_intersect_transform_traits<BaseTraits, Kernel>& traits)
    :m_traits(traits)
    {}

    template<typename Query>
    bool operator()(const Query& q, const Bounding_box& bbox) const
    {
      // TODO use base_traits
      return CGAL::do_intersect(q, compute_transformed_bbox(bbox));
    }

    template<typename Query>
    bool operator()(const Query& q, const Primitive& pr) const
    {
      internal_AABB::Transformed_datum_do_intersect<Kernel> f;
      return f(q, m_traits.transformation(), internal::Primitive_helper<BaseTraits>::get_datum(pr,m_traits));
    }

    // intersection with AABB-tree
    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Primitive& pr) const
    {
      return other_tree.do_intersect( internal::Primitive_helper<BaseTraits>::get_datum(pr,m_traits).transform(m_traits.transformation()));
    }

    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Bounding_box& bbox) const
    {
      return other_tree.do_intersect(compute_transformed_bbox(bbox));
    }
  };

  Do_intersect do_intersect_object() const{
    return Do_intersect(*this);
  }

  //Specific
  void set_transformation(const Aff_transformation_3<Kernel>& trans) const
  {
    m_transfo = trans;
  }

  const Aff_transformation_3<Kernel>& transformation() const { return m_transfo; }

};
}//end CGAL
#endif //CGAL_AABB_AABB_do_intersect_transform_traits_H
