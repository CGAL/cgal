// Copyright (c) 2010 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H
#define CGAL_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <boost/next_prior.hpp>
#include <list>
#include <vector>
#include <map>

namespace CGAL {

  template <class K>
  class Triangle_3;

namespace internal{

template <class Kernel>
void intersection_coplanar_triangles_cutoff(
  const typename Kernel::Point_3& p,
  const typename Kernel::Point_3& q,
  const typename Kernel::Point_3& r,
  const Kernel& k,
  std::list<typename Kernel::Point_3>& inter_pts
)
{
  if ( inter_pts.empty() ) return;
  typedef typename std::list<typename Kernel::Point_3>::iterator Iterator;
  typename Kernel::Coplanar_orientation_3 orient=k.coplanar_orientation_3_object();
  typename Kernel::Construct_line_3 Line_3=k.construct_line_3_object();
  //orient(p,q,r,r) is POSITIVE
  std::map<const typename Kernel::Point_3*,Orientation> orientations;
  for (Iterator it=inter_pts.begin();it!=inter_pts.end();++it)
    orientations[ &(*it) ]=orient(p,q,r,*it);

  int pt_added=0;
  
  const typename Kernel::Point_3* prev = &(*boost::prior(inter_pts.end()));
  Iterator stop = inter_pts.size() > 2 ? inter_pts.end() : boost::prior(inter_pts.end());
  for (Iterator it=inter_pts.begin();it!=stop;++it)
  {
    const typename Kernel::Point_3& curr=*it;
    Orientation or_prev=orientations[prev],or_curr=orientations[&curr];
    if ( (or_prev==POSITIVE && or_curr==NEGATIVE) || (or_prev==NEGATIVE && or_curr==POSITIVE) )
    {
      typename Intersection_traits<Kernel, typename Kernel::Line_3,
                                     typename Kernel::Line_3>
        ::result_type 
        obj = intersection(Line_3(p,q),Line_3(*prev,curr),k);
      // assert "not empty"
      CGAL_kernel_assertion(obj);
      const typename Kernel::Point_3* inter=intersect_get<typename Kernel::Point_3>(obj);
      CGAL_kernel_assertion(inter!=NULL);
      prev=&(* inter_pts.insert(it,*inter) );
      orientations[prev]=COLLINEAR;
      ++pt_added;
    }
    prev=&(*it);    
  }
  
  CGAL_kernel_assertion(pt_added<3);
  Iterator it=inter_pts.begin();
  while(it!=inter_pts.end())
  {
    if (orientations[&(*it)]==NEGATIVE)
      inter_pts.erase(it++);
    else
      ++it;
  }
}
  
template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3>::result_type
intersection_coplanar_triangles(
  const typename K::Triangle_3& t1,
  const typename K::Triangle_3& t2,
  const K& k)
{
  const typename K::Point_3& p=t1.vertex(0),q=t1.vertex(1),r=t1.vertex(2);
  
  std::list<typename K::Point_3> inter_pts;
  inter_pts.push_back(t2.vertex(0));
  inter_pts.push_back(t2.vertex(1));
  inter_pts.push_back(t2.vertex(2));

  //intersect t2 with the three half planes which intersection defines t1
  intersection_coplanar_triangles_cutoff(p,q,r,k,inter_pts); //line pq
  intersection_coplanar_triangles_cutoff(q,r,p,k,inter_pts); //line qr
  intersection_coplanar_triangles_cutoff(r,p,q,k,inter_pts); //line rp
  
  switch ( inter_pts.size() ) {
    case 0:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
    case 1:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>( * inter_pts.begin() );
    case 2:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>( k.construct_segment_3_object()(*inter_pts.begin(),
                                                         *boost::next(inter_pts.begin())) );
    case 3:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>( k.construct_triangle_3_object()(*inter_pts.begin(),
                                                          *boost::next(inter_pts.begin()),
                                                          *boost::prior(inter_pts.end())) );
    default:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>( std::vector<typename K::Point_3>(inter_pts.begin(),inter_pts.end()) );
  }
}

template<typename K>
struct Triangle_Line_visitor {
  typedef typename Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3 >
  ::result_type result_type;

  result_type
  operator()(const typename K::Point_3& p, const typename K::Segment_3& s) const {
    if ( s.has_on(p) )
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(p);
    else
      return result_type();
  }

  result_type
  operator()(const typename K::Segment_3& s, const typename K::Point_3& p) const {
    return operator()(p,s);
  }

  result_type
  operator()(const typename K::Point_3& p1, const typename K::Point_3& p2) const {
    if (p1==p2)
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(p1);
    else
      return result_type();
  }

  result_type
  operator()(const typename K::Segment_3& s1, const typename K::Segment_3& s2) const {
    typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type 
      v = intersection_collinear_segments(s1,s2,K());

    if(v) {
      if(const typename K::Segment_3* s = intersect_get<typename K::Segment_3>(v))
        return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(*s);
      if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
        return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(*p);
    }

    return result_type();
  }
};
  
template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3>::result_type
intersection(
  const typename K::Triangle_3& t1,
  const typename K::Triangle_3& t2,
  const K& k)
{
  CGAL_precondition(!t1.is_degenerate() && !t2.is_degenerate());

  typename Intersection_traits<K, typename K::Plane_3, typename K::Plane_3>::result_type 
    v = internal::intersection(t1.supporting_plane(), t2.supporting_plane(), k);

  if(!v) {
    // empty plane plane intersection
    return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
  }

  if(intersect_get<typename K::Plane_3>(v)) {
    //coplanar triangles
    return intersection_coplanar_triangles(t1,t2,k);
  }

  if(const typename K::Line_3* line=intersect_get<typename K::Line_3>(v)) {
    //The supporting planes of the triangles intersect along a line.
    typedef
    typename Intersection_traits<K, typename K::Triangle_3, typename K::Line_3>::result_type
      Triangle_Line_Inter;
    
    Triangle_Line_Inter inter1 = intersection_coplanar(t1,*line,k);
    Triangle_Line_Inter inter2 = intersection_coplanar(t2,*line,k);
    if(!inter1 || !inter2) {
      // one of the intersection is empty
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
    }
    #if CGAL_INTERSECTION_VERSION < 2
    // apply the binary visitor manually
    Triangle_Line_visitor<K> vis;
    if(const typename K::Point_3* p1 = intersect_get<typename K::Point_3>(inter1)) {
      if(const typename K::Point_3* p2 = intersect_get<typename K::Point_3>(inter2)) {
        return vis(*p1, *p2);
      } else if(const typename K::Segment_3* s2 = intersect_get<typename K::Segment_3>(inter2)) {
        return vis(*p1, *s2);
      }
    } else if(const typename K::Segment_3* s1 = intersect_get<typename K::Segment_3>(inter1)) {
      if(const typename K::Point_3* p2 = intersect_get<typename K::Point_3>(inter2)) {
        return vis(*s1, *p2);
      } else if(const typename K::Segment_3* s2 = intersect_get<typename K::Segment_3>(inter2)) {
        return vis(*s1, *s2);
      }
    }
    #else
    return boost::apply_visitor(Triangle_Line_visitor<K>(), *inter1, *inter2);
    #endif
  }
  return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
}

}//namespace internal

CGAL_INTERSECTION_FUNCTION_SELF(Triangle_3, 3)

} // namespace CGAL




#endif //CGAL_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H
