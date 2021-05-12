// Copyright (c) 2013 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_POINT_INSIDE_POLYHEDRON_RAY_3_TRIANGLE_3_TRAVERSAL_TRAITS_H
#define CGAL_POINT_INSIDE_POLYHEDRON_RAY_3_TRIANGLE_3_TRAVERSAL_TRAITS_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>


#include <boost/logic/tribool.hpp>
#include <CGAL/tags.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Intersections_3/Ray_3_Triangle_3.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/internal/AABB_tree/AABB_node.h>

namespace CGAL {
namespace internal {

template<typename AABBTraits, class Kernel, class Helper, class Tag_ray_is_vertical=Tag_false>
class Ray_3_Triangle_3_traversal_traits
{
protected:
  //the status indicates whether the query point is strictly inside the polyhedron, and the number of intersected triangles if yes
  std::pair<boost::logic::tribool,std::size_t>& m_status;
  bool m_stop;
  const AABBTraits& m_aabb_traits;
  typedef typename AABBTraits::Primitive Primitive;
  typedef CGAL::AABB_node<AABBTraits> Node;
  Helper m_helper;

public:
  Ray_3_Triangle_3_traversal_traits(std::pair<boost::logic::tribool,std::size_t>& status,
                                    const AABBTraits& aabb_traits,
                                    const Helper& h)
    :m_status(status),m_stop(false), m_aabb_traits(aabb_traits), m_helper(h)
  {m_status.first=true;}

  bool go_further() const { return !m_stop; }

  template<class Query>
  void intersection(const Query& query, const Primitive& primitive)
  {
    Intersections::internal::r3t3_do_intersect_endpoint_position_visitor visitor;
    std::pair<bool,Intersections::internal::R3T3_intersection::type> res=
      Intersections::internal::do_intersect(m_helper.get_primitive_datum(primitive, m_aabb_traits), query,Kernel(),visitor);

    if (res.first){
      switch (res.second){
        case Intersections::internal::R3T3_intersection::CROSS_FACET:
          ++m_status.second;
        break;
        case Intersections::internal::R3T3_intersection::ENDPOINT_IN_TRIANGLE:
          m_status.first=false;
          m_stop=true;
        break;
        default:
          m_status.first=boost::logic::indeterminate;
          m_stop=true;
      }
    }
  }

  template<class Query>
  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_aabb_traits.do_intersect_object()(query, m_helper.get_node_bbox(node));
  }
};


//specialization for vertical ray
template<typename AABBTraits, class Kernel, class TraversalTraits>
class Ray_3_Triangle_3_traversal_traits<AABBTraits,Kernel,TraversalTraits,Tag_true>:
  public Ray_3_Triangle_3_traversal_traits<AABBTraits,Kernel,TraversalTraits,Tag_false>
{
  typedef Ray_3_Triangle_3_traversal_traits<AABBTraits,Kernel,TraversalTraits,Tag_false> Base;
  typedef typename Kernel::Point_3 Point;
  typedef typename Base::Primitive Primitive;
  typedef CGAL::AABB_node<AABBTraits> Node;

public:
  Ray_3_Triangle_3_traversal_traits(std::pair<boost::logic::tribool,std::size_t>& status, const AABBTraits& aabb_traits, const TraversalTraits& tt)
    :Base(status, aabb_traits, tt){}

  template <class Query>
  bool do_intersect(const Query& query, const Bbox_3& bbox) const
  {
    const Point& source=query.point(0);
    const Point& target=query.point(1);

    bool inc_z=target.z()>source.z();

    //the ray does not intersect the z-slab
    if ( ( inc_z && source.z()>bbox.zmax() )|| (!inc_z && source.z()<bbox.zmin()) ) return false;

    //the source is not in the x-slab
    if (source.x() > bbox.xmax() || source.x()<bbox.xmin()) return false;
    //check if the source is not in the y-slab
    return source.y() <= bbox.ymax() && source.y()>=bbox.ymin();
  }

  template <class Query>
  bool do_intersect(const Query& query, const Node& node) const
  {
    return do_intersect(query,this->m_helper.get_node_bbox(node));
  }

private:
  typename Kernel::Point_2 x_project(const typename Kernel::Point_3& p) const{
    return typename Kernel::Point_2(p.y(),p.z());
  }
  typename Kernel::Point_2 y_project(const typename Kernel::Point_3& p) const{
    return typename Kernel::Point_2(p.x(),p.z());
  }
  typename Kernel::Point_2 z_project(const typename Kernel::Point_3& p) const{
    return typename Kernel::Point_2(p.x(),p.y());
  }
public:
  template<class Query>
  void intersection(const Query& query, const Primitive& primitive)
  {
    typename Kernel::Triangle_3 t = this->m_helper.get_primitive_datum(primitive, this->m_aabb_traits);
    if ( !do_intersect(query,t.bbox()) ) return;

    typename Kernel::Point_2 p0=z_project(t[0]);
    typename Kernel::Point_2 p1=z_project(t[1]);
    typename Kernel::Point_2 p2=z_project(t[2]);
    int indices[3]={0,1,2}; //to track whether triangle points have been swapt
    typename Kernel::Point_2 q=z_project( query.source() );

    Orientation orient_2=orientation(p0,p1,p2);

    //check whether the face has a normal vector in the xy-plane
    if (orient_2==COLLINEAR){
      //in that case the projection of the triangle along the z-axis is a segment.
      const typename Kernel::Point_2& other_point = p0!=p1?p1:p2;
      if ( orientation(p0,other_point,q) != COLLINEAR ) return;///no intersection

      //check if the ray source is above or below the triangle and compare it
      //with the direction of the ray
      //TODO and if yes return
      //this is just an optimisation, the current code is valid

      this->m_status.first=boost::logic::indeterminate;
      this->m_stop=true;
      return;
    }

    //regular case
    if (orient_2==NEGATIVE){
      std::swap(p1,p2);
      std::swap(indices[1],indices[2]);
    }

    //check whether the ray intersect the supporting plane
    Orientation orient_3 = orientation(t[indices[0]],t[indices[1]],t[indices[2]],query.source());
    if ( orient_3!=COPLANAR &&
          (
            //indicates whether the ray is oriented toward the positive side of the plane
            ( POSITIVE == CGAL::sign( query.to_vector().z() )  )
              ==
            //indicates whether the source of the ray is in the positive side of the plane
            (orient_3==POSITIVE)
          )
    ) return; //no intersection

    //position against first segment
    switch( orientation(p0,p1,q) ){
      case COLLINEAR:
        this->m_status.first=boost::logic::indeterminate;
        this->m_stop=true;
      case NEGATIVE:
        return;
      default:
      {}
    }
    //position against second segment
    switch( orientation(p1,p2,q) ){
      case COLLINEAR:
        this->m_status.first=boost::logic::indeterminate;
        this->m_stop=true;
      case NEGATIVE:
        return;
      default:
      {}
    }
    //position against third segment
    switch( orientation(p2,p0,q) ){
      case COLLINEAR:
        this->m_status.first=boost::logic::indeterminate;
        this->m_stop=true;
      case NEGATIVE:
        return;
      default:
      {}
    }

    if (orient_3==COPLANAR){
      //the endpoint is inside the triangle
      this->m_status.first=false;
      this->m_stop=true;
    }
    else
      ++(this->m_status.second);
  }
};

//special case when ray query is from another Kernel K1 is the kernel compatible with the AABB-tree
template<typename AABBTraits, class K1, class K2, class Helper>
class K2_Ray_3_K1_Triangle_3_traversal_traits
{
  //the status indicates whether the query point is strictly inside the polyhedron, and the number of intersected triangles if yes
  std::pair<boost::logic::tribool,std::size_t>& m_status;
  bool m_stop;
  const AABBTraits& m_aabb_traits;
  typedef typename AABBTraits::Primitive Primitive;
  typedef CGAL::AABB_node<AABBTraits> Node;
  Helper m_helper;
  CGAL::Cartesian_converter<K1,K2> to_K2;

public:
  K2_Ray_3_K1_Triangle_3_traversal_traits(std::pair<boost::logic::tribool,std::size_t>& status,
                                          const AABBTraits& aabb_traits,
                                          const Helper& h)
    :m_status(status), m_stop(false), m_aabb_traits(aabb_traits), m_helper(h)
  {m_status.first=true;}

  bool go_further() const { return !m_stop; }

  template<class Query>
  void intersection(const Query& query, const Primitive& primitive)
  {
    Intersections::internal::r3t3_do_intersect_endpoint_position_visitor visitor;
    std::pair<bool,Intersections::internal::R3T3_intersection::type> res=
      Intersections::internal::do_intersect(to_K2(m_helper.get_primitive_datum(primitive, m_aabb_traits)),
                                            query, K2(), visitor);

    if (res.first){
      switch (res.second){
        case Intersections::internal::R3T3_intersection::CROSS_FACET:
          ++m_status.second;
        break;
        case Intersections::internal::R3T3_intersection::ENDPOINT_IN_TRIANGLE:
          m_status.first=false;
          m_stop=true;
        break;
        default:
          m_status.first=boost::logic::indeterminate;
          m_stop=true;
      }
    }
  }

  template<class Query>
  bool do_intersect(const Query& query, const Node& node) const
  {
    return CGAL::do_intersect(query, m_helper.get_node_bbox(node));
  }
};


}// namespace internal
}// namespace CGAL

#endif
