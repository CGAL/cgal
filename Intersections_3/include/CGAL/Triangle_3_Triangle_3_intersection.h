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

#include <boost/next_prior.hpp>
#include <list>
#include <vector>
#include <map>

namespace CGAL {

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
      Object obj= intersection(Line_3(p,q),Line_3(*prev,curr),k);
      const typename Kernel::Point_3* inter=object_cast<typename Kernel::Point_3>(&obj);
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
  
template <class Kernel>
Object
intersection_coplanar_triangles(
  const typename Kernel::Triangle_3& t1,
  const typename Kernel::Triangle_3& t2,
  const Kernel& k)
{
  const typename Kernel::Point_3& p=t1.vertex(0),q=t1.vertex(1),r=t1.vertex(2);
  
  std::list<typename Kernel::Point_3> inter_pts;
  inter_pts.push_back(t2.vertex(0));
  inter_pts.push_back(t2.vertex(1));
  inter_pts.push_back(t2.vertex(2));

  //intersect t2 with the three half planes which intersection defines t1
  intersection_coplanar_triangles_cutoff(p,q,r,k,inter_pts); //line pq
  intersection_coplanar_triangles_cutoff(q,r,p,k,inter_pts); //line qr
  intersection_coplanar_triangles_cutoff(r,p,q,k,inter_pts); //line rp
  
  switch ( inter_pts.size() ) {
    case 0:
      return Object();
    case 1:
      return make_object( * inter_pts.begin() );
    case 2:
      return make_object( k.construct_segment_3_object()(*inter_pts.begin(),
                                                         *boost::next(inter_pts.begin())) );
    case 3:
      return make_object( k.construct_triangle_3_object()(*inter_pts.begin(),
                                                          *boost::next(inter_pts.begin()),
                                                          *boost::prior(inter_pts.end())) );
    default:
      return make_object( std::vector<typename Kernel::Point_3>(inter_pts.begin(),inter_pts.end()) );
  }
}
  
template <class Kernel>
Object
intersection(
  const typename Kernel::Triangle_3& t1,
  const typename Kernel::Triangle_3& t2,
  const Kernel& k)
{
  CGAL_precondition(!t1.is_degenerate() && !t2.is_degenerate());
  
  typename Kernel::Intersect_3 inter=k.intersect_3_object();
  CGAL::Object res=inter(t1.supporting_plane(),t2.supporting_plane());
  
  const typename Kernel::Line_3* line=CGAL::object_cast<typename Kernel::Line_3>(&res);
  
  if (line==NULL){
    const typename Kernel::Plane_3* plane=CGAL::object_cast<typename Kernel::Plane_3>(&res);
    if (plane!=NULL)
      return intersection_coplanar_triangles(t1,t2,k);
    return Object();
  }
  
  //The supporting planes of the triangles intersect along a line.
  Object inter1=intersection_coplanar(t1,*line,k);
  Object inter2=intersection_coplanar(t2,*line,k);
  
  
  const typename Kernel::Segment_3* sgt1=CGAL::object_cast<typename Kernel::Segment_3>(&inter1);
  
  if (sgt1 == NULL){
    //intersection of the line and triangle 1 is a point or is empty
    const typename Kernel::Point_3* pt1=CGAL::object_cast<typename Kernel::Point_3>(&inter1);
    if (pt1==NULL) return Object(); //the line does not intersect the triangle t1
    const typename Kernel::Segment_3* sgt2=CGAL::object_cast<typename Kernel::Segment_3>(&inter2);
    if (sgt2==NULL){
      const typename Kernel::Point_3* pt2=CGAL::object_cast<typename Kernel::Point_3>(&inter2);
      if (pt2==NULL || *pt1!=*pt2) return Object();
      return inter1;
    }
    //case point, segment
    if ( sgt2->has_on(*pt1) ) return inter1;
    return Object();
  }
  
  const typename Kernel::Segment_3* sgt2=CGAL::object_cast<typename Kernel::Segment_3>(&inter2);
  if (sgt2==NULL){
    const typename Kernel::Point_3* pt2=CGAL::object_cast<typename Kernel::Point_3>(&inter2);
    if ( pt2==NULL || !sgt1->has_on(*pt2) ) return Object();
    return inter2;
  }
  
  return intersection_collinear_segments(*sgt1,*sgt2,k);
}


}//namespace internal

template <class K>
inline
Object
intersection(const Triangle_3<K> &t1, const Triangle_3<K> &t2)
{
  return typename K::Intersect_3()(t1, t2);
}

} // namespace CGAL




#endif //CGAL_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H
