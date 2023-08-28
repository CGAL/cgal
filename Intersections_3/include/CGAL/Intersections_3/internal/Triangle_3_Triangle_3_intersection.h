// Copyright (c) 2010 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H

//#define CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Triangle_3_intersection.h>
#include <CGAL/Intersections_3/internal/Line_3_Line_3_intersection.h>
#include <CGAL/Intersections_3/internal/Plane_3_Plane_3_intersection.h>

#include <CGAL/kernel_assertions.h>

#include <boost/next_prior.hpp>
#include <boost/container/flat_set.hpp>

#include <list>
#include <map>
#include <vector>

namespace CGAL {
namespace Intersections {
namespace internal{

template <class Kernel>
struct Point_on_triangle
{
  // triangle points are not stored in this class but are expected
  // to always be passed in the same order. For a triangle pqr,
  // edge 0 is pq, edge 1 qr and edge 2 rp. Point 0 is p, 1 is q and 2 is r.
  //
  // (id, -1) point on t1
  // (-1, id) point on t2
  // (id1, id2) intersection of edges
  std::pair<int, int> t1_t2_ids;
  boost::container::flat_set<int> extra_t1; // store other ids of edges containing the point
  typename Kernel::FT alpha; //

//////

  static
  inline
  const typename Kernel::Point_3&
  point_from_id(const typename Kernel::Point_3& p,
                const typename Kernel::Point_3& q,
                const typename Kernel::Point_3& r,
                int id)
  {
    switch(id)
    {
      case 0:
        return p;
      case 1:
        return q;
      default:
        return r;
    }
  }

  Point_on_triangle(int i1, int i2=-1, typename Kernel::FT alpha = 0.) // TODO add global zero()?
    : t1_t2_ids(i1,i2)
    , alpha(alpha)
  {}

  // orientation of the current point wrt to edge id1 (p1q1)
  Orientation
  orientation (const typename Kernel::Point_3& p1, // source of edge edge_id1
               const typename Kernel::Point_3& q1, // target of edge edge_id1
               const typename Kernel::Point_3& r1,
               int edge_id1,
               const typename Kernel::Point_3& p2,
               const typename Kernel::Point_3& q2,
               const typename Kernel::Point_3& r2,
               const Kernel& k) const
  {
    if (t1_t2_ids.first!=-1)
    {
      if (t1_t2_ids.second==-1)
        return (edge_id1==t1_t2_ids.first || (edge_id1+1)%3==t1_t2_ids.first) ? ZERO:POSITIVE; // it is a point on t1
      // this is an intersection point

      if (t1_t2_ids.first==edge_id1)
        return ZERO;
      if (t1_t2_ids.first==(edge_id1+1)%3)
      {
        if (alpha==0) return ZERO;
        return alpha>=0 ? POSITIVE:NEGATIVE;
      }
      CGAL_assertion((t1_t2_ids.first+1)%3==edge_id1);
      if (alpha==1) return ZERO;
      return alpha<=1?POSITIVE:NEGATIVE;
    }
    else
    {
      //this is an input point of t2
      typename Kernel::Coplanar_orientation_3 orient = k.coplanar_orientation_3_object();
      const typename Kernel::Point_3& query = point_from_id(p2,q2,r2,t1_t2_ids.second);
      return orient(p1,q1,r1,query);
    }
  }

  int id1() const { return t1_t2_ids.first; }
  int id2() const { return t1_t2_ids.second; }

  // construct the intersection point from the info stored
  typename Kernel::Point_3
  point(const typename Kernel::Point_3& p1,
        const typename Kernel::Point_3& q1,
        const typename Kernel::Point_3& r1,
        const typename Kernel::Point_3& p2,
        const typename Kernel::Point_3& q2,
        const typename Kernel::Point_3& r2,
        const Kernel& k) const
  {
    if (t1_t2_ids.first==-1)
      return point_from_id(p2,q2,r2,t1_t2_ids.second);
    if (t1_t2_ids.second==-1)
      return point_from_id(p1,q1,r1,t1_t2_ids.first);

    return k.construct_barycenter_3_object()(point_from_id(p1,q1,r1,(t1_t2_ids.first+1)%3), alpha, point_from_id(p1,q1,r1,t1_t2_ids.first)) ;
  }
};

// the intersection of two triangles is computed by interatively intersection t2
// with halfspaces defined by edges of t1. The following function is called
// for each each on t1 on edge of the current intersection.
// pq is such an edge and p1q1 from t1 defines the halfspace intersection
// we are currently interseted in. We return the intersection point of
// pq with p1q1
template <class Kernel>
Point_on_triangle<Kernel>
intersection(const Point_on_triangle<Kernel>& p,
             const Point_on_triangle<Kernel>& q,
             int edge_id_t1,
             const typename Kernel::Point_3& p1,
             const typename Kernel::Point_3& q1,
//             const typename Kernel::Point_3& r1,
             const typename Kernel::Point_3& p2,
             const typename Kernel::Point_3& q2,
             const typename Kernel::Point_3& r2,
             const Kernel& k)
{
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "    calling intersection: ";
  std::cout << " (" << p.id1() << "," << p.id2() << ",[" << p.alpha << "]) -";
  std::cout << " (" << q.id1() << "," << q.id2() << ",[" << q.alpha << "]) || e" << edge_id_t1;
#endif
  typename Kernel::Compute_alpha_for_coplanar_triangle_intersection_3 compute_alpha
    = k.compute_alpha_for_coplanar_triangle_intersection_3_object();
  typedef Point_on_triangle<Kernel> Pot;
  switch(p.id1())
  {
    case -1:
    {
      switch(q.id1())
      {
        case -1: // A: (-1, ip2) - (-1, iq2)
        {
          CGAL_assertion((p.id2()+1)%3 == q.id2() || (q.id2()+1)%3 == p.id2());
//          CGAL_assertion(p.extra_t1.empty() && q.extra_t1.empty());  // TMP to see if it's worth implementing special case
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
          std::cout << " -- case 1\n";
#endif
          typename Kernel::FT alpha = compute_alpha(p1, q1,
                                                    Pot::point_from_id(p2, q2, r2, p.id2()),
                                                    Pot::point_from_id(p2, q2, r2, q.id2()));
          int id2 = (p.id2()+1)%3 == q.id2() ? p.id2() : q.id2();
          return  Point_on_triangle<Kernel>(edge_id_t1, id2, alpha); // intersection with an original edge of t2
        }
        default:
          if (q.id2()!=-1) // B: (-1, ip2) - (iq1, iq2)
          {
            if (p.id2() == q.id2() || p.id2() == (q.id2()+1)%3)
            {
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
              std::cout << " -- case 2\n";
#endif
              // points are on the same edge of t2 --> we shorten an already cut edge
              typename Kernel::FT alpha = compute_alpha(p1, q1,
                                                        Pot::point_from_id(p2, q2, r2,  q.id2()),
                                                        Pot::point_from_id(p2, q2, r2, (q.id2()+1)%3));

              return  Point_on_triangle<Kernel>(edge_id_t1, q.id2(), alpha);
            }
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
            std::cout << " -- case 3\n";
#endif
            // point of t1: look for an edge of t1 containing both points
            CGAL_assertion( p.extra_t1.count(q.id1())!=0 || p.extra_t1.count(3-q.id1()-edge_id_t1)!=0 );
            int eid1 = p.extra_t1.count(q.id1())!=0 ? q.id1() : 3-q.id1()-edge_id_t1;
            return  Point_on_triangle<Kernel>((eid1+1)%3==edge_id_t1?edge_id_t1:(edge_id_t1+1)%3, -1); // vertex of t1
          }
          // C: (-1, ip2) - (iq1, -1)
          //vertex of t1, special case t1 edge passed thru a vertex of t2
          CGAL_assertion(edge_id_t1 == 2);
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
          std::cout << " -- case 4\n";
#endif
          CGAL_assertion(q.id1()==1);
          CGAL_assertion(!p.extra_t1.empty());
          return  Point_on_triangle<Kernel>(p.extra_t1.count(0)==1?0:2,-1);
      }
    }
    default:
    {
      switch(p.id2())
      {
        case -1:
        {
          switch(q.id1())
          {
            case -1: // G: (ip1, -1) - (-1, iq2)
              //vertex of t1, special case t1 edge passed thru a vertex of t2
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
              std::cout << " -- case 5\n";
#endif
              CGAL_assertion(edge_id_t1 == 2);
              CGAL_assertion(p.id1()==1);
              CGAL_assertion(!q.extra_t1.empty());
              return  Point_on_triangle<Kernel>(q.extra_t1.count(0)==1?0:2,-1);
            default:
            {
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
              std::cout << " -- case 6\n";
#endif
              CGAL_assertion(q.id2()!=-1); // I: (ip1, -1) - (iq2, -1)
              //H: (ip1,-1), (iq1, iq2)
              CGAL_assertion(edge_id_t1==2);
              // p and q are on the same edge of t1
              CGAL_assertion(p.id1()==q.id1() || p.id1()==(q.id1()+1)%3);
              return  Point_on_triangle<Kernel>((q.id1()+1)%3==edge_id_t1?edge_id_t1:(edge_id_t1+1)%3 , -1);
            }
          }
        }
        default:
        {
          switch(q.id1())
          {
            case -1: // D: (ip1, ip2) - (-1, iq2)
            {
              if (q.id2() == p.id2() || q.id2() == (p.id2()+1)%3)
              {
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
                std::cout << " -- case 7\n";
#endif
                // points are on the same edge of t2 --> we shorten an already cut edge
                typename Kernel::FT alpha = compute_alpha(p1, q1,
                                                          Pot::point_from_id(p2, q2, r2,  p.id2()),
                                                          Pot::point_from_id(p2, q2, r2, (p.id2()+1)%3));

                return  Point_on_triangle<Kernel>(edge_id_t1, p.id2(), alpha);
              }
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
              std::cout << " -- case 8\n";
#endif
              // point of t1
              //std::cout << "q.extra_t1: "; for(int qet1 : q.extra_t1) std::cout << " " << qet1; std::cout << "\n";
              CGAL_assertion( q.extra_t1.count(p.id1())!=0 || q.extra_t1.count(3-p.id1()-edge_id_t1)!=0 );
              int eid1 = q.extra_t1.count(p.id1())!=0 ? p.id1() : 3-p.id1()-edge_id_t1;
              return  Point_on_triangle<Kernel>((eid1+1)%3==edge_id_t1?edge_id_t1:(edge_id_t1+1)%3, -1); // vertex of t1
            }
            default:
            {
              switch(q.id2())
              {
                case -1: // F: (ip1, ip2) - (iq1, -1)
                {
                  // p and q are on the same edge of t1
                  CGAL_assertion(q.id1()==p.id1() ||  q.id1()==(p.id1()+1)%3);
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
                  std::cout << " -- case 9\n";
#endif
                  return  Point_on_triangle<Kernel>((p.id1()+1)%3==edge_id_t1?edge_id_t1:(edge_id_t1+1)%3 , -1);
                }
                default: // E: (ip1, ip2) - (iq1, iq2)
                {
                  if (p.id2()==q.id2())
                  {
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
                    std::cout << " -- case 10\n";
#endif
                    typename Kernel::FT alpha = compute_alpha(p1, q1,
                                                              Pot::point_from_id(p2, q2, r2,  q.id2()),
                                                              Pot::point_from_id(p2, q2, r2, (q.id2()+1)%3));
                    return  Point_on_triangle<Kernel>(edge_id_t1, q.id2(), alpha);
                  }
                  // we are intersecting an edge of t1
                  CGAL_assertion(p.id1()==q.id1() || edge_id_t1==2);
                  int eid1 = p.id1()==q.id1() ? p.id1() : 1;
                  return  Point_on_triangle<Kernel>((eid1+1)%3==edge_id_t1?edge_id_t1:(edge_id_t1+1)%3, -1); // vertex of t1
                }
              }
            }
          }
        }
      }
    }
  }
}

template <class Kernel>
void intersection_coplanar_triangles_cutoff(const typename Kernel::Point_3& p1,
                                            const typename Kernel::Point_3& q1,
                                            const typename Kernel::Point_3& r1,
                                            int edge_id,
                                            const typename Kernel::Point_3& p2,
                                            const typename Kernel::Point_3& q2,
                                            const typename Kernel::Point_3& r2,
                                            const Kernel& k,
                                            std::list<Point_on_triangle<Kernel>>& inter_pts)
{
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  cutoff using e" << edge_id << ": "
            << to_double(p1.x()) << " " << to_double(p1.y()) << " " << to_double(p1.z()) << " "
            << to_double(q1.x()) << " " << to_double(q1.y()) << " " << to_double(q1.z()) << "\n";
#endif
  typedef typename std::list<Point_on_triangle<Kernel>>::iterator Iterator;

  if(inter_pts.empty())
    return;

  //orient(p1,q1,r1,r1) is POSITIVE
  std::map<const Point_on_triangle<Kernel>*,Orientation> orientations; // TODO skip map
  for (Point_on_triangle<Kernel>& pot : inter_pts)
  {
    orientations[ &pot ]=pot.orientation(p1,q1,r1,edge_id,p2,q2,r2,k);
    if (pot.id1()==-1 && orientations[ &pot ]==COLLINEAR)
      pot.extra_t1.insert(edge_id);
  }

#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "    Orientations:";
  for (const Point_on_triangle<Kernel>& pot : inter_pts)
    std::cout << " " << orientations[ &pot ];
  std::cout << "\n";
#endif
  CGAL_kernel_assertion_code(int pt_added = 0);

  Iterator prev = std::prev(inter_pts.end());

  Iterator stop = inter_pts.size() > 2 ? inter_pts.end() : std::prev(inter_pts.end());
  for(Iterator it=inter_pts.begin(); it!=stop; ++it)
  {
    Orientation or_prev = orientations[&(*prev)],
                or_curr = orientations[&(*it)];

    if((or_prev == POSITIVE && or_curr == NEGATIVE) ||
       (or_prev == NEGATIVE && or_curr == POSITIVE))
    {
      Point_on_triangle<Kernel> new_pt = intersection(*prev, *it, edge_id, p1, q1, p2, q2, r2, k);

      prev = inter_pts.insert(it,new_pt);
      orientations[&(*prev)] = COLLINEAR;
      CGAL_kernel_assertion_code(++pt_added);
    }

    prev = it;
  }

  CGAL_kernel_assertion(pt_added<3);
  Iterator it = inter_pts.begin();
  while(it!=inter_pts.end())
  {
    if(orientations[&(*it)] == NEGATIVE)
      inter_pts.erase(it++);
    else
      ++it;
  }
}

template <class K>
typename Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3>::result_type
intersection_coplanar_triangles(const typename K::Triangle_3& t1,
                                const typename K::Triangle_3& t2,
                                const K& k)
{
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  auto to_string = [](const typename K::Triangle_3& t)
  {
    std::stringstream sstr;
    sstr << "4 "
         << to_double(t[0].x()) << " "  << to_double(t[0].y()) << " "  << to_double(t[0].z()) << " "
         << to_double(t[1].x()) << " "  << to_double(t[1].y()) << " "  << to_double(t[1].z()) << " "
         << to_double(t[2].x()) << " "  << to_double(t[2].y()) << " "  << to_double(t[2].z()) << " "
         << to_double(t[0].x()) << " "  << to_double(t[0].y()) << " "  << to_double(t[0].z()) << "\n";
    return sstr.str();
  };

  std::cout << "intersection_coplanar_triangles\n";
  std::ofstream("/tmp/t1.polylines.txt") << std::setprecision(17) << to_string(t1) << "\n";
  std::ofstream("/tmp/t2.polylines.txt") << std::setprecision(17) << to_string(t2) << "\n";
#endif
  const typename K::Point_3& p1 = t1.vertex(0),
                             q1 = t1.vertex(1),
                             r1 = t1.vertex(2);

  const typename K::Point_3& p2 = t2.vertex(0),
                             q2 = t2.vertex(1),
                             r2 = t2.vertex(2);

  std::list<Point_on_triangle<K>> inter_pts;
  inter_pts.push_back(Point_on_triangle<K>(-1,0));
  inter_pts.push_back(Point_on_triangle<K>(-1,1));
  inter_pts.push_back(Point_on_triangle<K>(-1,2));

#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  auto print_points = [&]()
  {
    for(auto p : inter_pts) std::cout << "  (" << p.id1() << "," << p.id2() << ",[" << p.alpha << "]) "; std::cout <<"\n";
  };
  std::cout << "  ipts size: " << inter_pts.size() << "\n";
  print_points();
#endif
  //intersect t2 with the three half planes which intersection defines t1
  intersection_coplanar_triangles_cutoff(p1,q1,r1,0,p2,q2,r2,k,inter_pts); //line pq
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << inter_pts.size() << "\n";
  print_points();
#endif
  intersection_coplanar_triangles_cutoff(q1,r1,p1,1,p2,q2,r2,k,inter_pts); //line qr
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << inter_pts.size() << "\n";
  print_points();
#endif
  intersection_coplanar_triangles_cutoff(r1,p1,q1,2,p2,q2,r2,k,inter_pts); //line rp
#ifdef CGAL_DEBUG_COPLANAR_T3_T3_INTERSECTION
  std::cout << "  ipts size: " << inter_pts.size() << "\n";
  print_points();
#endif

  auto point = [&](const Point_on_triangle<K>& pot){ return pot.point(p1,q1,r1,p2,q2,r2,k); };
  switch(inter_pts.size())
  {
    case 0:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
    case 1:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(point(*inter_pts.begin()));
    case 2:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(
            k.construct_segment_3_object()(point(*inter_pts.begin()), point(*std::next(inter_pts.begin()))) );
    case 3:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(
            k.construct_triangle_3_object()(point(*inter_pts.begin()), point(*std::next(inter_pts.begin())), point(*std::prev(inter_pts.end()))) );
    default:
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(
            std::vector<typename K::Point_3>(boost::make_transform_iterator(inter_pts.begin(), point),
                                             boost::make_transform_iterator(inter_pts.end(), point)));
  }
}

template<typename K>
struct Triangle_Line_visitor
{
  typedef typename Intersection_traits<K, typename K::Triangle_3, typename K::Triangle_3 >::result_type result_type;

  result_type operator()(const typename K::Point_3& p, const typename K::Segment_3& s) const
  {
    if(typename K::Has_on_3()(s, p))
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(p);
    else
      return result_type();
  }

  result_type operator()(const typename K::Segment_3& s, const typename K::Point_3& p) const
  {
    return operator()(p,s);
  }

  result_type operator()(const typename K::Point_3& p1, const typename K::Point_3& p2) const
  {
    if(p1 == p2)
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>(p1);
    else
      return result_type();
  }

  result_type operator()(const typename K::Segment_3& s1, const typename K::Segment_3& s2) const
  {
    typename Intersection_traits<K, typename K::Segment_3, typename K::Segment_3>::result_type
      v = intersection_collinear_segments(s1,s2,K());

    if(v)
    {
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
intersection(const typename K::Triangle_3& t1,
             const typename K::Triangle_3& t2,
             const K& k)
{
  CGAL_precondition(!t1.is_degenerate() && !t2.is_degenerate());

  typename Intersection_traits<K, typename K::Plane_3, typename K::Plane_3>::result_type
    v = internal::intersection(t1.supporting_plane(), t2.supporting_plane(), k);

  if(!v)
  {
    // empty plane plane intersection
    return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
  }

  if(intersect_get<typename K::Plane_3>(v))
  {
    //coplanar triangles
    return intersection_coplanar_triangles(t1,t2,k);
  }

  if(const typename K::Line_3* line=intersect_get<typename K::Line_3>(v))
  {
    //The supporting planes of the triangles intersect along a line.
    typedef typename Intersection_traits<K, typename K::Triangle_3, typename K::Line_3>::result_type Triangle_Line_Inter;

    Triangle_Line_Inter inter1 = intersection_coplanar(t1,*line,k);
    Triangle_Line_Inter inter2 = intersection_coplanar(t2,*line,k);
    if(!inter1 || !inter2)
    {
      // one of the intersection is empty
      return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
    }

    return std::visit(Triangle_Line_visitor<K>(), *inter1, *inter2);
  }

  return intersection_return<typename K::Intersect_3, typename K::Triangle_3, typename K::Triangle_3>();
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_TRIANGLE_3_TRIANGLE_3_INTERSECTION_H
