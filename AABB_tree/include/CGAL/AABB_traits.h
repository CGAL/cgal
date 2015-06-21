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
// $URL$
// $Id$
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_TRAITS_H_
#define CGAL_AABB_TRAITS_H_

#include <CGAL/Bbox_3.h>
#include <CGAL/AABB_intersections.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
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

  #if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
  template <typename ... T>
  void set_shared_data(T&& ... t){
    m_primitive_data=Primitive::construct_shared_data(std::forward<T>(t)...);
  }
  #else
  void set_shared_data(){
    m_primitive_data=Primitive::construct_shared_data();
  }

  template <class T1>
  void set_shared_data(T1& t1){
    m_primitive_data=Primitive::construct_shared_data(t1);
  }

  template <class T1,class T2>
  void set_shared_data(T1& t1, T2& t2){
    m_primitive_data=Primitive::construct_shared_data(t1,t2);
  }

  template <class T1,class T2,class T3>
  void set_shared_data(T1& t1,T2& t2,T3& t3){
    m_primitive_data=Primitive::construct_shared_data(t1,t2,t3);
  }

  template <class T1,class T2,class T3,class T4>
  void set_shared_data(T1& t1,T2& t2,T3& t3,T4& t4){
    m_primitive_data=Primitive::construct_shared_data(t1,t2,t3,t4);
  }

  template <class T1,class T2,class T3,class T4,class T5>
  void set_shared_data(T1& t1,T2& t2,T3& t3,T4& t4,T5& t5){
    m_primitive_data=Primitive::construct_shared_data(t1,t2,t3,t4,t5);
  }
  #endif
  const typename Primitive::Shared_data& shared_data() const {return m_primitive_data;}
};

} } //end of namespace internal::AABB_tree

/// \addtogroup PkgAABB_tree
/// @{

/// This traits class handles any type of 3D geometric
/// primitives provided that the proper intersection tests and
/// constructions are implemented. It handles points, rays, lines and
/// segments as query types for intersection detection and
/// computations, and it handles points as query type for distance
/// queries.
/// \cgalModels AABBTraits
/// \tparam GeomTraits must  be a model of the concept \ref AABBGeomTraits,
/// snd provide the geometric types as well as the intersection tests and computations.
/// \tparam Primitive provide the type of primitives stored in the AABB_tree.
///   It is a model of the concept `AABBPrimitive` or `AABBPrimitiveWithSharedData`.
///
/// \sa `AABBTraits`
/// \sa `AABB_tree`
/// \sa `AABBPrimitive`
/// \sa `AABBPrimitiveWithSharedData`
template<typename GeomTraits, typename AABBPrimitive>
class AABB_traits:
  public internal::AABB_tree::AABB_traits_base<AABBPrimitive>
{
  typedef typename CGAL::Object Object;
public:
  typedef AABB_traits<GeomTraits, AABBPrimitive> AT;
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

  ///
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


  /// Default constructor.
  AABB_traits() { };


  typedef typename GeomTraits::Compute_squared_distance_3 Squared_distance;
  Squared_distance squared_distance_object() const { return GeomTraits().compute_squared_distance_3_object(); }

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
  class Sort_primitives
  {
    const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
  public:
    Sort_primitives(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
      : m_traits(traits) {}

    template<typename PrimitiveIterator>
    void operator()(PrimitiveIterator first,
                    PrimitiveIterator beyond,
                    const typename AT::Bounding_box& bbox) const
      {
        PrimitiveIterator middle = first + (beyond - first)/2;
        switch(longest_axis(bbox))
        {
        case AT::CGAL_AXIS_X: // sort along x
          std::nth_element(first, middle, beyond, boost::bind(less_x,_1,_2,m_traits));
          break;
        case AT::CGAL_AXIS_Y: // sort along y
          std::nth_element(first, middle, beyond, boost::bind(less_y,_1,_2,m_traits));
          break;
        case AT::CGAL_AXIS_Z: // sort along z
          std::nth_element(first, middle, beyond, boost::bind(less_z,_1,_2,m_traits));
          break;
        default:
          CGAL_error();
        }
      }
  };

  Sort_primitives sort_primitives_object() const {return Sort_primitives(*this);}


  /*
   * Computes the bounding box of a set of primitives
   * @param first an iterator on the first primitive
   * @param beyond an iterator on the past-the-end primitive
   * @return the bounding box of the primitives of the iterator range
   */
  class Compute_bbox {
    const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
  public:
    Compute_bbox(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
      :m_traits (traits) {}

    template<typename ConstPrimitiveIterator>
    typename AT::Bounding_box operator()(ConstPrimitiveIterator first,
                                         ConstPrimitiveIterator beyond) const
      {
        typename AT::Bounding_box bbox = compute_bbox(*first,m_traits);
        for(++first; first != beyond; ++first)
        {
          bbox = bbox + compute_bbox(*first,m_traits);
        }
        return bbox;
      }
  };

  Compute_bbox compute_bbox_object() const {return Compute_bbox(*this);}


  class Do_intersect {
    const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
  public:
    Do_intersect(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
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
  };

  Do_intersect do_intersect_object() const {return Do_intersect(*this);}

class Intersection {
  const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
public:
  Intersection(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
    :m_traits(traits) {}
    #if CGAL_INTERSECTION_VERSION < 2
template<typename Query>
boost::optional<typename AT::Object_and_primitive_id>
operator()(const Query& query, const typename AT::Primitive& primitive) const
{
  typedef boost::optional<Object_and_primitive_id> Intersection;

  CGAL::Object object = GeomTraits().intersect_3_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
  if ( object.empty() )
    return Intersection();
  else
    return Intersection(Object_and_primitive_id(object,primitive.id()));
}
    #else
    template<typename Query>
    boost::optional< typename Intersection_and_primitive_id<Query>::Type >
    operator()(const Query& query, const typename AT::Primitive& primitive) const {
      typename cpp11::result_of<typename GeomTraits::Intersect_3(Query, typename Primitive::Datum) >::type
        inter_res = GeomTraits().intersect_3_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
      if (!inter_res)
          return boost::optional<typename Intersection_and_primitive_id<Query>::Type>();
      return boost::make_optional( std::make_pair(*inter_res, primitive.id()) );
    }
    #endif
};

Intersection intersection_object() const {return Intersection(*this);}

  // This should go down to the GeomTraits, i.e. the kernel
  class Closest_point {
      typedef typename AT::Point_3 Point;
      typedef typename AT::Primitive Primitive;
    const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
  public:
    Closest_point(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
      : m_traits(traits) {}


    Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
    {
        return CGAL::nearest_point_3(p, internal::Primitive_helper<AT>::get_datum(pr,m_traits), bound);
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
  static Bounding_box compute_bbox (const Primitive& pr,
                                    const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  {
    return internal::Primitive_helper<AT>::get_datum(pr,traits).bbox();
  }

  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  static Axis longest_axis(const Bounding_box& bbox);

  /// Comparison functions
  static bool less_x(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).x() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).x(); }
  static bool less_y(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).y() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).y(); }
  static bool less_z(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).z() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).z(); }

};  // end class AABB_traits


//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
template<typename GT, typename P>
typename AABB_traits<GT,P>::Axis
AABB_traits<GT,P>::longest_axis(const Bounding_box& bbox)
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

#endif // CGAL_AABB_TRAITS_H_
