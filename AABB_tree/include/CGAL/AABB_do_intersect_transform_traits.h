
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

namespace CGAL {
// forward declaration
template< typename AABBTraits>
class AABB_tree;

namespace internal_AABB
{
template<class Kernel>
struct Actual_intersect
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
  struct Filter_test
  {
    template<class Query, class Datum_t>
    bool operator()(const Query& q,
                    const CGAL::Aff_transformation_3<Kernel>& transfo,
                    const Datum_t& pr) const
    {
      return Actual_intersect<Kernel>()(q,transfo,pr);
    }
  };
  template<class Kernel>
  struct Filter_test<Kernel, Tag_true>
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
      
      typedef internal_AABB::Actual_intersect<EK> Exactator;
      typedef internal_AABB::Actual_intersect<FK> Approxator;
      typedef CGAL::Filtered_predicate<
          Exactator, 
          Approxator, 
          C2E, 
          C2F> Filter;
      Filter fi;
      return fi(q, transfo,pr);
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
    typedef Simple_cartesian<Interval_nt_advanced>             Approximate_kernel;
    typedef Cartesian_converter<Kernel, Approximate_kernel>    C2F;
    
    
    const AABB_do_intersect_transform_traits<BaseTraits, Kernel>& m_traits;
    C2F c2f;
  public:
    Do_intersect(const AABB_do_intersect_transform_traits<BaseTraits, Kernel>& traits)
    :m_traits(traits)
    {}
    
    template<typename Query>
    bool operator()(const Query& q, const Bounding_box& bbox) const
    {
      Point_3 min(bbox.xmin(), bbox.ymin(), bbox.zmin()),
          max(bbox.xmax(), bbox.ymax(), bbox.zmax());
      
      typename Approximate_kernel::Point_3 app_min, 
          app_max;
        Bounding_box temp_box =
        min.bbox() + max.bbox();
        Point_3 tmin(temp_box.xmin(), temp_box.ymin(), temp_box.zmin()),
            tmax(temp_box.xmax(), temp_box.ymax(), temp_box.zmax());
        app_min=c2f(tmin);
        app_max = c2f(tmax);
      Bounding_box transfo_box(to_double(app_min.x().inf()), to_double(app_min.y().inf()), to_double(app_min.z().inf()),
                               to_double(app_max.x().sup()), to_double(app_max.y().sup()), to_double(app_max.z().sup()));
      
      bool res = CGAL::do_intersect(c2f(q), transfo_box);
      return res;
    }
    
    template<typename Query>
    bool operator()(const Query& q, const Primitive& pr) const
    {
      internal_AABB::Filter_test<Kernel> f;
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
      Point_3 min(bbox.xmin(), bbox.ymin(), bbox.zmin()),
          max(bbox.xmax(), bbox.ymax(), bbox.zmax());
      
      min = m_traits.transformation().transform(min);
      max = m_traits.transformation().transform(max);
      
      Bounding_box temp_box =
          min.bbox() + max.bbox();
      min=Point_3(temp_box.xmin(), temp_box.ymin(), temp_box.zmin());
      max=Point_3(temp_box.xmax(), temp_box.ymax(), temp_box.zmax());
      Bounding_box transfo_box(to_double(min.x()), to_double(min.y()), to_double(min.z()),
                               to_double(max.x()), to_double(max.y()), to_double(max.z()));
      return other_tree.do_intersect(transfo_box);
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
