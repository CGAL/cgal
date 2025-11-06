// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_TRAITS_2_H_
#define CGAL_AABB_TRAITS_2_H_

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Default.h>
#include <CGAL/intersections.h>
#include <CGAL/AABB_tree/internal/AABB_traits_base.h>
#include <CGAL/AABB_tree/internal/Has_nested_type_Shared_data.h>
#include <CGAL/AABB_tree/internal/Is_ray_intersection_geomtraits.h>
#include <CGAL/AABB_tree/internal/Primitive_helper.h>
#include <CGAL/AABB_tree/internal/Remove_optional.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <CGAL/Kernel_23/internal/Has_boolean_tags.h>
#include <CGAL/Search_traits_2.h>
#include <optional>

/// \file AABB_traits_2.h

namespace CGAL {

namespace internal{  namespace AABB_tree {


// AABB_traits_intersection_base_2 brings in the Intersection_distance predicate,
// if GeomTraits is a model RayIntersectionGeomTraits.
template <typename GeomTraits, bool ray_intersection_geom_traits=Is_ray_intersection_geomtraits_2<GeomTraits>::value>
struct AABB_traits_intersection_base_2;

template <typename GeomTraits>
struct AABB_traits_intersection_base_2<GeomTraits,false>{};

template <typename GeomTraits>
struct AABB_traits_intersection_base_2<GeomTraits,true>{
  template<typename AABBTree, typename SkipFunctor>
  friend class AABB_ray_intersection;
private:
  typedef typename GeomTraits::FT    FT;
  typedef typename GeomTraits::Point_2 Point;
  typedef typename GeomTraits::Cartesian_const_iterator_2 Cartesian_const_iterator;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator;

  static Construct_cartesian_const_iterator construct_cartesian_const_iterator_object() {
    return GeomTraits().construct_cartesian_const_iterator_2_object();
  }

public:
  typedef typename GeomTraits::Ray_2 Ray;
  typedef typename GeomTraits::Vector_2 Vector;
  typedef typename GeomTraits::Construct_source_2 Construct_source;
  typedef typename GeomTraits::Construct_vector_2 Construct_vector;

  static Construct_source construct_source_object() {
    return GeomTraits().construct_source_2_object();
  }

  static Construct_vector construct_vector_object() {
    return GeomTraits().construct_vector_2_object();
  }

  // Defining Bounding_box and other types from the full AABB_traits_2
  // here might seem strange, but otherwise we would need to use
  // CRTP to get access to the derived class, which would bloat the
  // code more.
  typedef typename CGAL::Bbox_2      Bounding_box;

  struct Intersection_distance {
    std::optional<FT> operator()(const Ray& ray, const Bounding_box& bbox) const {
      FT t_near = -DBL_MAX; // std::numeric_limits<FT>::lowest(); C++1903
      FT t_far = DBL_MAX;

      const Construct_cartesian_const_iterator construct_cartesian_const_iterator_2
        = GeomTraits().construct_cartesian_const_iterator_2_object();
      const Construct_source construct_source_2 = GeomTraits().construct_source_2_object();
      const Construct_vector construct_vector_2 = GeomTraits().construct_vector_2_object();
      const Point source = construct_source_2(ray);
      const Vector direction = construct_vector_2(ray);
      Cartesian_const_iterator source_iter = construct_cartesian_const_iterator_2(source);
      Cartesian_const_iterator direction_iter = construct_cartesian_const_iterator_2(direction);

      for(int i = 0; i < 2; ++i, ++source_iter, ++direction_iter) {
        if(*direction_iter == 0) {
          if((*source_iter < (bbox.min)(i)) || (*source_iter > (bbox.max)(i))) {
            return std::nullopt;
          }
        } else {
          FT t1 = ((bbox.min)(i) - *source_iter) / *direction_iter;
          FT t2 = ((bbox.max)(i) - *source_iter) / *direction_iter;

          t_near = (std::max)(t_near, (std::min)(t1, t2));
          t_far = (std::min)(t_far, (std::max)(t1, t2));

          if(t_near > t_far || t_far < FT(0.))
            return std::nullopt;
        }
      }

      if(t_near < FT(0.))
        return FT(0.);
      else
        return t_near;
    }
  };

  Intersection_distance intersection_distance_object() const { return Intersection_distance(); }
};

template<typename GeomTraits>
class Compare_distance_2 {
  typedef typename GeomTraits::Point_2 Point;
  typedef typename GeomTraits::FT FT;

  /// Bounding box type.
  typedef typename CGAL::Bbox_2 Bounding_box;
public:
  CGAL::Comparison_result operator()(const Point& p, const Bounding_box& bb, const Point& bound) const
  {
    return do_intersect_circle_iso_rectangle_2
    (GeomTraits().construct_circle_2_object()
      (p, GeomTraits().compute_squared_distance_2_object()(p, bound)), bb) ?
      CGAL::SMALLER : CGAL::LARGER;
  }

  template <class Solid>
  CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const Point& bound) const
  {
    return GeomTraits().do_intersect_2_object()
      (GeomTraits().construct_circle_2_object()
        (p, GeomTraits().compute_squared_distance_2_object()(p, bound)), pr) ?
      CGAL::SMALLER : CGAL::LARGER;
  }

  template <class Solid>
  CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const FT& sq_distance) const
  {
    return GeomTraits().do_intersect_2_object()
      (GeomTraits().construct_circle_2_object()(p, sq_distance),
        pr) ?
      CGAL::SMALLER :
      CGAL::LARGER;
  }

  typename GeomTraits::Boolean do_intersect_circle_iso_rectangle_2(const typename GeomTraits::Circle_2& circle,
    const typename GeomTraits::Iso_rectangle_2& rec) const
  {
    typedef typename GeomTraits::FT       FT;
    typedef typename GeomTraits::Point_2  Point;

    Point center = circle.center();

    // Check that the minimum distance to the box is smaller than the radius, otherwise there is
    // no intersection. `distance` stays at 0 if the center is inside or on `rec`.
    FT distance = FT(0);
    if (center.x() < rec.xmin())
    {
      FT d = rec.xmin() - center.x();
      distance += d * d;
    }
    else if (center.x() > rec.xmax())
    {
      FT d = center.x() - rec.xmax();
      distance += d * d;
    }

    if (center.y() < rec.ymin())
    {
      FT d = rec.ymin() - center.y();
      distance += d * d;
    }
    else if (center.y() > rec.ymax())
    {
      FT d = center.y() - rec.ymax();
      distance += d * d;
    }

    if (distance <= circle.squared_radius())
      return true;

    return false;
  }
};

template <typename GeomTraits, bool Has_filtered_predicates = internal::Has_filtered_predicates<GeomTraits>::value, bool Has_static_filters = internal::Has_static_filters<GeomTraits>::value>
class Compare_distance_getter_2 {};

template <typename GeomTraits>
class Compare_distance_getter_2<GeomTraits, false, false> {
  // this class is in charge of checking what K provides (i.e., can we use filtered predicates, can we use statically filtered predicates, etc.)
  // depending on that it defines
public:
  typedef Compare_distance_2<GeomTraits> type;
  static Compare_distance_2<GeomTraits> compare_distance_object() {
    return Compare_distance_2<GeomTraits>();
  }
};

template <typename GeomTraits>
class Compare_distance_getter_2<GeomTraits, true, false> {
  // this class is in charge of checking what K provides (i.e., can we use filtered predicates, can we use statically filtered predicates, etc.)
  // depending on that it defines

  typedef GeomTraits                                              Kernel;

  typedef typename Kernel::Exact_kernel                           EKernel;
  typedef typename Kernel::Approximate_kernel                     AKernel;
  typedef typename Kernel::C2E                                    C2E;
  typedef typename Kernel::C2F                                    C2F;

  typedef Compare_distance_2<EKernel> Exact_functor;
  typedef Compare_distance_2<AKernel> Filtered_functor;

public:
  typedef Filtered_predicate<Exact_functor, Filtered_functor,
    C2E, C2F> Compare_distance_pred;
  typedef Compare_distance_pred type;

  static Compare_distance_pred compare_distance_object() {
    return Compare_distance_pred(Exact_functor(), Filtered_functor());
  }
};

template <typename GeomTraits>
class Compare_distance_getter_2<GeomTraits, true, true> {
  // this class is in charge of checking what K provides (i.e., can we use filtered predicates, can we use statically filtered predicates, etc.)
  // depending on that it defines
  class Statically_filtered_compare_distance {
  public:
    typedef typename GeomTraits::Point_2 Point;
    typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Circle_2 Circle_2;

    /// Bounding box type.
    typedef CGAL::Bbox_2 Bounding_box;

    template <class Solid>
    CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const Point& bound) const {
      return Compare_distance_getter_2<GeomTraits, true, false>::compare_distance_object()(p, pr, bound);
    }

    template <class Solid>
    CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const FT& sq_distance) const {
      return Compare_distance_getter_2<GeomTraits, true, false>::compare_distance_object()(p, pr, sq_distance);
    }

    Comparison_result operator()(const Point& p, const Bounding_box& b, const Point& bound) const {
      Circle_2 s = GeomTraits().construct_circle_2_object()(p, GeomTraits().compute_squared_distance_2_object()(p, bound));

      CGAL_BRANCH_PROFILER_3(std::string("semi-static failures/attempts/calls to   : ") +
        std::string(CGAL_PRETTY_FUNCTION), tmp);

      internal::Static_filters_predicates::Get_approx<Point> get_approx; // Identity functor for all points
      const Point& c = s.center();

      double scx, scy, ssr;
      double bxmin = b.xmin(), bymin = b.ymin(),
        bxmax = b.xmax(), bymax = b.ymax();

      if (internal::fit_in_double(get_approx(c).x(), scx) &&
        internal::fit_in_double(get_approx(c).y(), scy) &&
        internal::fit_in_double(s.squared_radius(), ssr))
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

        if ((ssr < 1.11261183279326254436e-293) || (ssr > 2.80889552322236673473e+306)) {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Compare_distance_getter_2<GeomTraits, true, false>::compare_distance_object()(p, b, bound);
        }
        double distance = 0;
        double max1 = 0;
        double double_tmp_result = 0;
        double eps = 0;
        if (scx < bxmin)
        {
          double bxmin_scx = bxmin - scx;
          max1 = bxmin_scx;

          distance = square(bxmin_scx);
          double_tmp_result = (distance - ssr);

          if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }
        else if (scx > bxmax)
        {
          double scx_bxmax = scx - bxmax;
          max1 = scx_bxmax;

          distance = square(scx_bxmax);
          double_tmp_result = (distance - ssr);

          if ((max1 < 3.33558365626356687717e-147) || (max1 > 1.67597599124282407923e+153))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }

        if (scy < bymin)
        {
          double bymin_scy = bymin - scy;
          if (max1 < bymin_scy) {
            max1 = bymin_scy;
          }

          distance += square(bymin_scy);
          double_tmp_result = (distance - ssr);

          if ((max1 < 3.33558365626356687717e-147) || ((max1 > 1.67597599124282407923e+153)))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps) {
            return CGAL::LARGER;
          }
        }
        else if (scy > bymax)
        {
          double scy_bymax = scy - bymax;
          if (max1 < scy_bymax) {
            max1 = scy_bymax;
          }
          distance += square(scy_bymax);
          double_tmp_result = (distance - ssr);

          if (((max1 < 3.33558365626356687717e-147)) || ((max1 > 1.67597599124282407923e+153)))
            return CGAL::SMALLER;

          eps = 1.99986535548615598560e-15 * (std::max)(ssr, square(max1));

          if (double_tmp_result > eps)
            return CGAL::LARGER;
        }

        // double_tmp_result and eps were growing all the time
        // no need to test for > eps as done earlier in at least one case

        return CGAL::SMALLER;
      }
      return Compare_distance_getter_2<GeomTraits, true, false>::compare_distance_object()(p, b, bound);
    }
  };
public:
  typedef Statically_filtered_compare_distance type;

  static Statically_filtered_compare_distance compare_distance_object() {
    return Statically_filtered_compare_distance();
  }
};

} } //end of namespace internal::AABB_tree

/// \addtogroup PkgAABBTreeRef
/// @{

// forward declaration
template< typename AABBTraits>
class AABB_tree;

template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default>
class AABB_traits_2
#ifndef DOXYGEN_RUNNING
: public internal::AABB_tree::AABB_traits_base<AABBPrimitive>,
  public internal::AABB_tree::AABB_traits_intersection_base_2<GeomTraits>,
  public Search_traits_2<GeomTraits>
#endif
{
  typedef typename CGAL::Object Object;
  typedef GeomTraits Geom_traits;
public:

  typedef AABB_traits_2<GeomTraits, AABBPrimitive, BboxMap> AT;
  // AABBTraits concept types
  typedef typename GeomTraits::FT FT;
  typedef AABBPrimitive Primitive;

  typedef typename std::pair<Object,typename Primitive::Id> Object_and_primitive_id;

  typedef typename std::pair<typename GeomTraits::Point_2, typename Primitive::Id> Point_and_primitive_id;

  /// `Intersection_and_primitive_id<Query>::%Type::first_type` is found according to
  /// the result type of `GeomTraits::Intersect_2::operator()`. If it is
  /// `std::optional<T>` then it is `T`, and the result type otherwise.
  template<typename Query>
  struct Intersection_and_primitive_id {
    typedef decltype(
      std::declval<typename GeomTraits::Intersect_2>()(
        std::declval<Query>(),
        std::declval<typename Primitive::Datum>())) Intersection_type;

    typedef std::pair<
      typename internal::Remove_optional<Intersection_type>::type,
      typename Primitive::Id > Type;
  };

  // types for search tree
  /// \name Types
  /// @{

  /// <summary>
  /// point type
  /// </summary>
  typedef typename GeomTraits::Point_2 Point;

  /// additional types for the search tree, required by the RangeSearchTraits concept
  /// \bug This is not documented for now in the AABBTraits concept.
  typedef typename GeomTraits::Iso_rectangle_2 Iso_rectangle_2;

  /// Bounding box type.
  typedef typename CGAL::Bbox_2 Bounding_box;

  /// @}

  typedef typename GeomTraits::Circle_2 Circle_2;
  typedef typename GeomTraits::Cartesian_const_iterator_2 Cartesian_const_iterator_2;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator_2;
  typedef typename GeomTraits::Construct_center_2 Construct_center_2;
  typedef typename GeomTraits::Compute_squared_radius_2 Compute_squared_radius_2;
  typedef typename GeomTraits::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef typename GeomTraits::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename GeomTraits::Construct_iso_rectangle_2 Construct_iso_rectangle_2;

  BboxMap bbm;

  /// Default constructor.
  AABB_traits_2() { }

  AABB_traits_2(BboxMap bbm)
    : bbm(bbm)
  {}


  typedef typename GeomTraits::Compute_squared_distance_2 Squared_distance;
  Squared_distance squared_distance_object() const { return GeomTraits().compute_squared_distance_2_object(); }

  typedef typename GeomTraits::Equal_2 Equal;
  Equal equal_object() const { return GeomTraits().equal_2_object(); }

  /**
   * @internal
   * @brief Sorts [first,beyond[
   * @param first iterator on first element
   * @param beyond iterator on beyond element
   * @param bbox the bounding box of [first,beyond[
   *
   * Sorts the range defined by [first,beyond[. Sort is achieved on bbox longest
   * axis, using the comparison function `<dim>_less_than` (dim in {x,y,z})
   */
  class Split_primitives
  {
    typedef AABB_traits_2<GeomTraits,AABBPrimitive,BboxMap> Traits;
    const Traits& m_traits;
  public:
    Split_primitives(const AABB_traits_2<GeomTraits,AABBPrimitive,BboxMap>& traits)
      : m_traits(traits) {}

    typedef void result_type;
    template<typename PrimitiveIterator>
    void operator()(PrimitiveIterator first,
                    PrimitiveIterator beyond,
                    const typename AT::Bounding_box& bbox) const
      {
        PrimitiveIterator middle = first + (beyond - first)/2;
        switch(Traits::longest_axis(bbox))
        {
        case AT::CGAL_AXIS_X: // sort along x
          std::nth_element(first, middle, beyond, [this](const Primitive& p1, const Primitive& p2){ return Traits::less_x(p1, p2, this->m_traits); });
          break;
        case AT::CGAL_AXIS_Y: // sort along y
          std::nth_element(first, middle, beyond, [this](const Primitive& p1, const Primitive& p2){ return Traits::less_y(p1, p2, this->m_traits); });
          break;
        default:
          CGAL_error();
        }
      }
  };

  Split_primitives split_primitives_object() const {return Split_primitives(*this);}


  /*
   * Computes the bounding box of a set of primitives
   * @param first an iterator on the first primitive
   * @param beyond an iterator on the past-the-end primitive
   * @return the bounding box of the primitives of the iterator range
   */
  class Compute_bbox {
    const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& m_traits;
  public:
    Compute_bbox(const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& traits)
      :m_traits (traits) {}

    template<typename ConstPrimitiveIterator>
    typename AT::Bounding_box operator()(ConstPrimitiveIterator first,
                                         ConstPrimitiveIterator beyond) const
    {
      typename AT::Bounding_box bbox = m_traits.compute_bbox(*first,m_traits.bbm);
      for(++first; first != beyond; ++first)
        {
          bbox = bbox + m_traits.compute_bbox(*first,m_traits.bbm);
        }
      return bbox;
    }

  };

  Compute_bbox compute_bbox_object() const {return Compute_bbox(*this);}

  /// \brief Function object using `GeomTraits::Do_intersect`.
  /// In the case the query is a `CGAL::AABB_tree`, the `do_intersect()`
  /// function of this tree is used.
  class Do_intersect {
    const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& m_traits;
  public:
    Do_intersect(const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& traits)
      :m_traits(traits) {}

    template<typename Query>
    bool operator()(const Query& q, const Bounding_box& bbox) const
    {
      return GeomTraits().do_intersect_2_object()(q, bbox);
    }

    template<typename Query>
    bool operator()(const Query& q, const Primitive& pr) const
    {
      return GeomTraits().do_intersect_2_object()(q, internal::Primitive_helper<AT>::get_datum(pr,m_traits));
    }

    // intersection with AABB-tree
    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Primitive& pr) const
    {
      return other_tree.do_intersect( internal::Primitive_helper<AT>::get_datum(pr,m_traits) );
    }

    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Bounding_box& bbox) const
    {
      return other_tree.do_intersect(bbox);
    }
  };

  Do_intersect do_intersect_object() const {return Do_intersect(*this);}


  class Intersection {
    const AABB_traits_2<GeomTraits,AABBPrimitive,BboxMap>& m_traits;
  public:
    Intersection(const AABB_traits_2<GeomTraits,AABBPrimitive,BboxMap>& traits)
      :m_traits(traits) {}
    template<typename Query>
    std::optional< typename Intersection_and_primitive_id<Query>::Type >
    operator()(const Query& query, const typename AT::Primitive& primitive) const {
      auto inter_res = GeomTraits().intersect_2_object()(query, internal::Primitive_helper<AT>::get_datum(primitive,m_traits));
      if (!inter_res)
        return std::nullopt;
      return std::make_optional( std::make_pair(*inter_res, primitive.id()) );
    }
  };

  Intersection intersection_object() const {return Intersection(*this);}


  // This should go down to the GeomTraits, i.e. the kernel
  class Closest_point {
      typedef typename AT::Point Point;
      typedef typename AT::Primitive Primitive;
    const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& m_traits;
  public:
    Closest_point(const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& traits)
      : m_traits(traits) {}


    Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
    {
      GeomTraits geom_traits;
      Point closest_point = geom_traits.construct_projected_point_2_object()(
        internal::Primitive_helper<AT>::get_datum(pr,m_traits), p);

      return (geom_traits.compare_distance_2_object()(p, closest_point, bound) == LARGER) ?
               bound : closest_point;
    }
  };

  typedef typename internal::AABB_tree::Compare_distance_getter_2<GeomTraits>::type Compare_distance;

  Closest_point closest_point_object() const { return Closest_point(*this); }
  Compare_distance compare_distance_object() const { return internal::AABB_tree::Compare_distance_getter_2<GeomTraits>::compare_distance_object(); }

  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1} Axis;

  static Axis longest_axis(const Bounding_box& bbox);

private:
  /**
   * @brief Computes bounding box of one primitive
   * @param pr the primitive
   * @return the bounding box of the primitive \c pr
   */
  template <typename PM>
  Bounding_box compute_bbox(const Primitive& pr, const PM&)const
  {
    return get(bbm, pr.id());
  }

  Bounding_box compute_bbox(const Primitive& pr, const Default&)const
  {
    return GeomTraits().construct_bbox_2_object()(internal::Primitive_helper<AT>::get_datum(pr, *this));
  }

  /// Comparison functions
  static bool less_x(const Primitive& pr1, const Primitive& pr2,const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& traits)
  {
    return GeomTraits().less_x_2_object()( internal::Primitive_helper<AT>::get_reference_point(pr1,traits),
                                           internal::Primitive_helper<AT>::get_reference_point(pr2,traits) );
  }
  static bool less_y(const Primitive& pr1, const Primitive& pr2,const AABB_traits_2<GeomTraits,AABBPrimitive, BboxMap>& traits)
  {
    return GeomTraits().less_y_2_object()( internal::Primitive_helper<AT>::get_reference_point(pr1,traits),
                                           internal::Primitive_helper<AT>::get_reference_point(pr2,traits) );
  }

};  // end class AABB_traits_2


//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
  template<typename GT, typename P, typename B>
  typename AABB_traits_2<GT,P,B>::Axis
    AABB_traits_2<GT,P,B>::longest_axis(const Bounding_box& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();

  if(dx>=dy)
    {
      return CGAL_AXIS_X;
    }
    else
    {
      return CGAL_AXIS_Y;
    }
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_TRAITS_2_H_
