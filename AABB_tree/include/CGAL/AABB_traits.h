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
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_TRAITS_H_
#define CGAL_AABB_TRAITS_H_

#include <CGAL/Bbox_3.h>
#include <CGAL/AABB_intersections.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <boost/optional.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/bind.hpp>

namespace CGAL {

  
namespace internal{

//for backward compatibility (if auto is available, use it)
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
template <class Primitive, bool has_shared_data=Has_nested_type_Shared_data<Primitive>::value>
struct Primitive_helper;
  
template <class Primitive>
struct Primitive_helper<Primitive,false>{
  typename Datum_result_type<Primitive>::type get_datum(const Primitive& p) const {return p.datum();}
  typename Point_result_type<Primitive>::type get_reference_point(const Primitive& p) const {return p.reference_point();}
};

template <class Primitive>
struct Primitive_helper<Primitive,true>{
  typename  Primitive::Shared_data m_primitive_data;
  
  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <class PrimitiveType, typename ... T>
  void set_shared_data(T ... t){
    m_primitive_data=PrimitiveType::construct_shared_data(t...);
  }  
  #else
  template <class PrimitiveType>
  void set_shared_data(){
    m_primitive_data=PrimitiveType::construct_shared_data();
  }

  template <class PrimitiveType, class T1>
  void set_shared_data(T1 t1){
    m_primitive_data=PrimitiveType::construct_shared_data(t1);
  }
  
  template <class PrimitiveType, class T1,class T2,class T3>
  void set_shared_data(T1 t1,T2 t2,T3 t3){
    m_primitive_data=PrimitiveType::construct_shared_data(t1,t2,t3);
  }

  template <class PrimitiveType, class T1,class T2,class T3,class T4>
  void set_shared_data(T1 t1,T2 t2,T3 t3,T4 t4){
    m_primitive_data=PrimitiveType::construct_shared_data(t1,t2,t3,t4);
  }
  
  template <class PrimitiveType, class T1,class T2,class T3,class T4,class T5>
  void set_shared_data(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5){
    m_primitive_data=PrimitiveType::construct_shared_data(t1,t2,t3,t4,t5);
  }
  #endif
  
  typename Datum_result_type<Primitive>::type get_datum(const Primitive& p) const {return p.datum(m_primitive_data);}
  typename Point_result_type<Primitive>::type get_reference_point(const Primitive& p) const {return p.reference_point(m_primitive_data);}
};

}
  
/**
 * @class AABB_traits
 *
 *
 */
template<typename GeomTraits, typename AABBPrimitive>
class AABB_traits:
  public internal::Primitive_helper<AABBPrimitive>
{
public:
  typedef AABB_traits<GeomTraits, AABBPrimitive> AT;
  /// AABBTraits concept types
  typedef typename CGAL::Bbox_3 Bounding_box;
  typedef typename CGAL::Object Object;

  typedef AABBPrimitive Primitive;
  typedef typename AABBPrimitive::Datum Datum;

  typedef typename GeomTraits::Point_3 Point;

  typedef typename std::pair<Object,typename Primitive::Id> Object_and_primitive_id;
  typedef typename std::pair<Point,typename Primitive::Id> Point_and_primitive_id;

  // types for search tree
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Sphere_3 Sphere_3;
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename GeomTraits::Construct_center_3 Construct_center_3;
  typedef typename GeomTraits::Construct_iso_cuboid_3 Construct_iso_cuboid_3;
  typedef typename GeomTraits::Construct_min_vertex_3 Construct_min_vertex_3;
  typedef typename GeomTraits::Construct_max_vertex_3 Construct_max_vertex_3;
  typedef typename GeomTraits::Compute_squared_radius_3 Compute_squared_radius_3;
  typedef typename GeomTraits::Compute_squared_distance_3 Compute_squared_distance_3;
  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_3
                     Construct_cartesian_const_iterator_3;

  /// Constructor
  AABB_traits() { };

  /// Non-virtual Destructor
  ~AABB_traits() { };


  /// 
  /**
   * @brief Sorts [first,beyond[
   * @param first iterator on first element
   * @param beyond iterator on beyond element
   * @param bbox the bounding box of [first,beyond[
   *
   * Sorts the range defined by [first,beyond[. Sort is achieved on bbox longuest
   * axis, using the comparison function <dim>_less_than (dim in {x,y,z})
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


  /**
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
      return GeomTraits().do_intersect_3_object()(q, m_traits.get_datum(pr));
    }
  };

  Do_intersect do_intersect_object() const {return Do_intersect(*this);}

  class Intersection {
    const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
  public:
    Intersection(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
      :m_traits(traits) {}

    template<typename Query>
    boost::optional<typename AT::Object_and_primitive_id>
    operator()(const Query& query, const typename AT::Primitive& primitive) const
    {
      typedef boost::optional<Object_and_primitive_id> Intersection;

      CGAL::Object object = GeomTraits().intersect_3_object()(m_traits.get_datum(primitive),query);
      if ( object.empty() )
        return Intersection();
      else
        return Intersection(Object_and_primitive_id(object,primitive.id()));
    }
  };

  Intersection intersection_object() const {return Intersection(*this);}


  // This should go down to the GeomTraits, i.e. the kernel
  class Closest_point {
    typedef typename AT::Point Point;
    typedef typename AT::Primitive Primitive;
    const AABB_traits<GeomTraits,AABBPrimitive>& m_traits;
  public:
    Closest_point(const AABB_traits<GeomTraits,AABBPrimitive>& traits)
      : m_traits(traits) {}


    Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
    {
        return CGAL::nearest_point_3(p, m_traits.get_datum(pr), bound);
    }
  };

  // This should go down to the GeomTraits, i.e. the kernel
  // and the internal implementation should change its name from
  // do_intersect to something like does_contain (this is what we compute,
  // this is not the same do_intersect as the spherical kernel)
  class Compare_distance {
      typedef typename AT::Point Point;
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
    return traits.get_datum(pr).bbox();
  }

  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  static Axis longest_axis(const Bounding_box& bbox);
  /// Comparison functions
  static bool less_x(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  { return traits.get_reference_point(pr1).x() < traits.get_reference_point(pr2).x(); }
  static bool less_y(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  { return traits.get_reference_point(pr1).y() < traits.get_reference_point(pr2).y(); }
  static bool less_z(const Primitive& pr1, const Primitive& pr2,const AABB_traits<GeomTraits,AABBPrimitive>& traits)
  { return traits.get_reference_point(pr1).z() < traits.get_reference_point(pr2).z(); }

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


}  // end namespace CGAL

#endif // CGAL_AABB_TRAITS_H_
