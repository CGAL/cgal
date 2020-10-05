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

#ifndef CGAL_AABB_TRAITS_H_
#define CGAL_AABB_TRAITS_H_

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>
#include <CGAL/intersections.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>

#include <boost/optional.hpp>
#include <boost/bind.hpp>

/// \file AABB_traits.h

namespace CGAL {

namespace internal{  namespace AABB_tree {

template <class T>
struct Remove_optional  { typedef T type; };

template <class T>
struct Remove_optional< ::boost::optional<T> >  { typedef T type; };

//helper controlling whether extra data should be stored in the AABB_tree traits class
template <class Primitive, bool has_shared_data=Has_nested_type_Shared_data<Primitive>::value>
struct AABB_traits_base;

template <class Primitive>
struct AABB_traits_base<Primitive,false>{};

template <class Primitive>
struct AABB_traits_base<Primitive,true>{
  typename  Primitive::Shared_data m_primitive_data;

  template <typename ... T>
  void set_shared_data(T&& ... t){
    m_primitive_data=Primitive::construct_shared_data(std::forward<T>(t)...);
  }
  const typename Primitive::Shared_data& shared_data() const {return m_primitive_data;}
};

// AABB_traits_base_2 brings in the Intersection_distance predicate,
// if GeomTraits is a model RayIntersectionGeomTraits.
template <typename GeomTraits, bool ray_intersection_geom_traits=Is_ray_intersection_geomtraits<GeomTraits>::value>
struct AABB_traits_base_2;

template <typename GeomTraits>
struct AABB_traits_base_2<GeomTraits,false>{};

template <typename GeomTraits>
struct AABB_traits_base_2<GeomTraits,true>{
  typedef typename GeomTraits::Ray_3 Ray_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;
  typedef typename GeomTraits::FT    FT;
  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_source_3 Construct_source_3;
  typedef typename GeomTraits::Construct_vector_3 Construct_vector_3;

  // Defining Bounding_box and other types from the full AABB_traits
  // here is might seem strange, but otherwise we would need to use
  // CRTP to get access to the derived class, which would bloat the
  // code more.
  typedef typename CGAL::Bbox_3      Bounding_box;

  struct Intersection_distance {
    boost::optional<FT> operator()(const Ray_3& ray, const Bounding_box& bbox) const {
      FT t_near = -DBL_MAX; // std::numeric_limits<FT>::lowest(); C++1903
      FT t_far = DBL_MAX;

      const Construct_cartesian_const_iterator_3 construct_cartesian_const_iterator_3
        = GeomTraits().construct_cartesian_const_iterator_3_object();
      const Construct_source_3 construct_source_3 = GeomTraits().construct_source_3_object();
      const Construct_vector_3 construct_vector_3 = GeomTraits().construct_vector_3_object();
      const Point_3 source = construct_source_3(ray);
      const Vector_3 direction = construct_vector_3(ray);
      Cartesian_const_iterator_3 source_iter = construct_cartesian_const_iterator_3(source);
      Cartesian_const_iterator_3 direction_iter = construct_cartesian_const_iterator_3(direction);

      for(int i = 0; i < 3; ++i, ++source_iter, ++direction_iter) {
        if(*direction_iter == 0) {
          if((*source_iter < (bbox.min)(i)) || (*source_iter > (bbox.max)(i))) {
            return boost::none;
          }
        } else {
          FT t1 = ((bbox.min)(i) - *source_iter) / *direction_iter;
          FT t2 = ((bbox.max)(i) - *source_iter) / *direction_iter;

          t_near = (std::max)(t_near, (std::min)(t1, t2));
          t_far = (std::min)(t_far, (std::max)(t1, t2));

          // if(t1 > t2)
          //   std::swap(t1, t2);
          // if(t1 > t_near)
          //   t_near = t1;
          // if(t2 < t_far)
          //   t_far = t2;

          if(t_near > t_far || t_far < FT(0.))
            return boost::none;
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

} } //end of namespace internal::AABB_tree

/// \addtogroup PkgAABBTreeRef
/// @{

// forward declaration
template< typename AABBTraits>
class AABB_tree;


/// This traits class handles any type of 3D geometric
/// primitives provided that the proper intersection tests and
/// constructions are implemented. It handles points, rays, lines and
/// segments as query types for intersection detection and
/// computations, and it handles points as query type for distance
/// queries.
///
/// \cgalModels AABBTraits
/// \cgalModels AABBRayIntersectionTraits

/// \tparam GeomTraits must  be a model of the concept \ref AABBGeomTraits,
/// and provide the geometric types as well as the intersection tests and computations.
/// \tparam Primitive provide the type of primitives stored in the AABB_tree.
///   It is a model of the concept `AABBPrimitive` or `AABBPrimitiveWithSharedData`.
///
/// \tparam BboxMap must be a model of `ReadablePropertyMap` that has as key type a primitive id,
///                 and as value type a `Bounding_box`.
///                 If the type is `Default` the `Datum` must have the
///                 member function `bbox()` that returns the bounding box  of the primitive.
///
/// If the argument `GeomTraits` is a model of the concept \ref
/// AABBRayIntersectionGeomTraits, this class is also a model of \ref
/// AABBRayIntersectionTraits.
///
/// \sa `AABBTraits`
/// \sa `AABB_tree`
/// \sa `AABBPrimitive`
/// \sa `AABBPrimitiveWithSharedData`

  template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default>
class AABB_traits
#ifndef DOXYGEN_RUNNING
: public internal::AABB_tree::AABB_traits_base<AABBPrimitive>,
  public internal::AABB_tree::AABB_traits_base_2<GeomTraits>
#endif
{
  typedef typename CGAL::Object Object;
public:
  typedef GeomTraits Geom_traits;

  typedef AABB_traits<GeomTraits, AABBPrimitive, BboxMap> AT;
  // AABBTraits concept types
  typedef typename GeomTraits::FT FT;
  typedef AABBPrimitive Primitive;

  typedef typename std::pair<Object,typename Primitive::Id> Object_and_primitive_id;

  typedef typename std::pair<typename GeomTraits::Point_3, typename Primitive::Id> Point_and_primitive_id;

  /// `Intersection_and_primitive_id<Query>::%Type::first_type` is found according to
  /// the result type of `GeomTraits::Intersect_3::operator()`,
  /// (that is cpp11::result_of<GeomTraits::Intersect_3(Query, Primitive::Datum)>::type). If it is
  /// `boost::optional<T>` then it is `T`, and the result type otherwise.
  template<typename Query>
  struct Intersection_and_primitive_id {
    typedef typename cpp11::result_of<
      typename GeomTraits::Intersect_3(Query, typename Primitive::Datum)
    >::type Intersection_type;

    typedef std::pair<
      typename internal::AABB_tree::Remove_optional<Intersection_type>::type,
      typename Primitive::Id > Type;
  };

  // types for search tree
  /// \name Types
  /// @{

  /// Point query type.
  typedef typename GeomTraits::Point_3 Point_3;

  /// additionnal types for the search tree, required by the RangeSearchTraits concept
  /// \bug This is not documented for now in the AABBTraits concept.
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

  /// Bounding box type.
  typedef typename CGAL::Bbox_3 Bounding_box;

  /// @}

  typedef typename GeomTraits::Sphere_3 Sphere_3;
  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_center_3 Construct_center_3;
  typedef typename GeomTraits::Compute_squared_radius_3 Compute_squared_radius_3;
  typedef typename GeomTraits::Construct_min_vertex_3 Construct_min_vertex_3;
  typedef typename GeomTraits::Construct_max_vertex_3 Construct_max_vertex_3;
  typedef typename GeomTraits::Construct_iso_cuboid_3 Construct_iso_cuboid_3;

  BboxMap bbm;

  /// Default constructor.
  AABB_traits() { }

  AABB_traits(BboxMap bbm)
    : bbm(bbm)
  {}


  typedef typename GeomTraits::Compute_squared_distance_3 Squared_distance;
  Squared_distance squared_distance_object() const { return GeomTraits().compute_squared_distance_3_object(); }

  typedef typename GeomTraits::Equal_3 Equal_3;
  Equal_3 equal_3_object() const { return GeomTraits().equal_3_object(); }

  /**
   * @internal
   * @brief Sorts [first,beyond[
   * @param first iterator on first element
   * @param beyond iterator on beyond element
   * @param bbox the bounding box of [first,beyond[
   *
   * Sorts the range defined by [first,beyond[. Sort is achieved on bbox longuest
   * axis, using the comparison function `<dim>_less_than` (dim in {x,y,z})
   */
  class Split_primitives
  {
    typedef AABB_traits<GeomTraits,AABBPrimitive,BboxMap> Traits;
    const Traits& m_traits;
  public:
    Split_primitives(const AABB_traits<GeomTraits,AABBPrimitive,BboxMap>& traits)
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
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_x,_1,_2,m_traits));
          break;
        case AT::CGAL_AXIS_Y: // sort along y
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_y,_1,_2,m_traits));
          break;
        case AT::CGAL_AXIS_Z: // sort along z
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_z,_1,_2,m_traits));
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
    const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& m_traits;
  public:
    Compute_bbox(const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& traits)
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
    const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& m_traits;
  public:
    Do_intersect(const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& traits)
      :m_traits(traits) {}

    template<typename Query>
    bool operator()(const Query& q, const Bounding_box& bbox) const
    {
      return CGAL::do_intersect(q, bbox);
    }

    template<typename Query>
    bool operator()(const Query& q, const Primitive& pr) const
    {
      return GeomTraits().do_intersect_3_object()(q, internal::Primitive_helper<AT>::get_datum(pr,m_traits));
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
    const AABB_traits<GeomTraits,AABBPrimitive,BboxMap>& m_traits;
  public:
    Intersection(const AABB_traits<GeomTraits,AABBPrimitive,BboxMap>& traits)
      :m_traits(traits) {}
    template<typename Query>
    boost::optional< typename Intersection_and_primitive_id<Query>::Type >
    operator()(const Query& query, const typename AT::Primitive& primitive) const {
      typename cpp11::result_of<typename GeomTraits::Intersect_3(Query, typename Primitive::Datum) >::type
        inter_res = GeomTraits().intersect_3_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
      if (!inter_res)
        return boost::none;
      return boost::make_optional( std::make_pair(*inter_res, primitive.id()) );
    }
  };

  Intersection intersection_object() const {return Intersection(*this);}


  // This should go down to the GeomTraits, i.e. the kernel
  class Closest_point {
      typedef typename AT::Point_3 Point;
      typedef typename AT::Primitive Primitive;
    const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& m_traits;
  public:
    Closest_point(const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& traits)
      : m_traits(traits) {}


    Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
    {
      GeomTraits geom_traits;
      Point closest_point = geom_traits.construct_projected_point_3_object()(
        internal::Primitive_helper<AT>::get_datum(pr,m_traits), p);
      return
        geom_traits.compare_distance_3_object()(p, closest_point, bound)==LARGER ?
        bound : closest_point;
    }
  };

  // This should go down to the GeomTraits, i.e. the kernel
  // and the internal implementation should change its name from
  // do_intersect to something like does_contain (this is what we compute,
  // this is not the same do_intersect as the spherical kernel)
  class Compare_distance {
      typedef typename AT::Point_3 Point;
      typedef typename AT::FT FT;
      typedef typename AT::Primitive Primitive;
  public:
      template <class Solid>
      CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const Point& bound) const
      {
          return GeomTraits().do_intersect_3_object()
          (GeomTraits().construct_sphere_3_object()
          (p, GeomTraits().compute_squared_distance_3_object()(p, bound)), pr)?
          CGAL::SMALLER : CGAL::LARGER;
      }

      template <class Solid>
      CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const FT& sq_distance) const
      {
        return GeomTraits().do_intersect_3_object()
          (GeomTraits().construct_sphere_3_object()(p, sq_distance),
           pr) ?
          CGAL::SMALLER :
          CGAL::LARGER;
      }
  };

  Closest_point closest_point_object() const {return Closest_point(*this);}
  Compare_distance compare_distance_object() const {return Compare_distance();}


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
    return internal::Primitive_helper<AT>::get_datum(pr,*this).bbox();
  }


  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  static Axis longest_axis(const Bounding_box& bbox);

  /// Comparison functions
  static bool less_x(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& traits)
  {
    return GeomTraits().less_x_3_object()( internal::Primitive_helper<AT>::get_reference_point(pr1,traits),
                                           internal::Primitive_helper<AT>::get_reference_point(pr2,traits) );
  }
  static bool less_y(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& traits)
  {
    return GeomTraits().less_y_3_object()( internal::Primitive_helper<AT>::get_reference_point(pr1,traits),
                                           internal::Primitive_helper<AT>::get_reference_point(pr2,traits) );
  }
  static bool less_z(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive, BboxMap>& traits)
  {
    return GeomTraits().less_z_3_object()( internal::Primitive_helper<AT>::get_reference_point(pr1,traits),
                                           internal::Primitive_helper<AT>::get_reference_point(pr2,traits) );
  }

};  // end class AABB_traits


//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
  template<typename GT, typename P, typename B>
  typename AABB_traits<GT,P,B>::Axis
  AABB_traits<GT,P,B>::longest_axis(const Bounding_box& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double dz = bbox.zmax() - bbox.zmin();

  if(dx>=dy)
  {
    if(dx>=dz)
    {
      return CGAL_AXIS_X;
    }
    else // dz>dx and dx>=dy
    {
      return CGAL_AXIS_Z;
    }
  }
  else // dy>dx
  {
    if(dy>=dz)
    {
      return CGAL_AXIS_Y;
    }
    else  // dz>dy and dy>dx
    {
      return CGAL_AXIS_Z;
    }
  }
}

/// @}

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_TRAITS_H_
