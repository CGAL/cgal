// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_INTERSECTION_DATA_STRUCTURE_H
#define CGAL_INTERSECTION_DATA_STRUCTURE_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Bbox_3.h>

#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <CGAL/intersections.h>

#include <numeric>

namespace CGAL {

template <class GT, class Type>
class Intersection_data_structure_3 {
  typedef Simple_cartesian<double> Double_kernel;
  typedef typename Double_kernel::Point_3 DPoint_3;

  typedef Segment_tree_map_traits_3<Double_kernel,
                                    Type> Traits;

  typedef typename GT::Triangle_3 Triangle_3;
  typedef typename GT::Segment_3 Segment_3;
  typedef typename GT::Point_3 Point_3;
public:
  typedef typename Traits::Interval Interval;
  typedef typename Traits::Pure_interval Pure_interval;
  typedef typename Traits::Key Key;

  enum Intersection_type { INSIDE, ON_BOUNDARY, NO_INTERSECTION };

  template <class T>
  Pure_interval make_pure_interval(const T& t) const
  {
    const Bbox_3& b = t.bbox();
    return std::make_pair(DPoint_3(b.xmin(), b.ymin(), b.zmin()),
                          DPoint_3(b.xmax()+epsilon,
                                   b.ymax()+epsilon,
                                   b.zmax()+epsilon));
  }

  Interval make_interval(const Type& t) const
  {
    return std::make_pair(make_pure_interval(t), t);
  }

  Object proper_intersection(const Triangle_3& t, const Segment_3& s) const
  {
    typename GT::Construct_supporting_plane_3 plane =
      gt.construct_supporting_plane_3_object();
    typename GT::Intersect_3 intersect =
      gt.intersect_3_object();


    if(do_intersect_properly(t, s) != INSIDE)
      return Object();
    return intersect(s, plane(t));
  }

  Object proper_intersection(const Segment_3& s, const Triangle_3& t) const
  {
    return proper_intersection(t, s);
  }

  Intersection_type do_intersect_properly(const Triangle_3& t, const Segment_3& s) const
  {
    typename GT::Construct_vertex_3 vertex =
      gt.construct_vertex_3_object();
    typename GT::Coplanar_3 coplanar =
      gt.coplanar_3_object();
    typename GT::Do_intersect_3 do_intersect =
      gt.do_intersect_3_object();

    const Point_3& p = vertex(s,0);
    const Point_3& q = vertex(s,1);

    const Point_3& a = vertex(t,0);
    const Point_3& b = vertex(t,1);
    const Point_3& c = vertex(t,2);

    if(coplanar(p, q, a,  b) ||
       coplanar(p, q, b,  c) ||
       coplanar(p, q, c,  a))
      return ON_BOUNDARY;
    return do_intersect(t, s) ? INSIDE : NO_INTERSECTION;
  }

  Intersection_type do_intersect_properly(const Segment_3& s, const Triangle_3& t) const
  {
    return do_intersect_properly(t, s);
  }

private:
  typedef Segment_tree_3<Traits> Segment_tree;

  typedef std::vector<Type> Elements;
  Elements elements;
  typename Elements::size_type nb_of_elements;
  Segment_tree tree;
  Bbox_3 bounding_box;
  double max_width;
  double epsilon;
  GT gt;

public:
  typedef Bbox_3 Bbox;

  Intersection_data_structure_3(GT gt = GT())
    : elements(),
      nb_of_elements(),
      tree(),
      bounding_box(0., 0., 0.,
                   0., 0., 0.),
      max_width(),
      epsilon(),
      gt(gt)
  {
  }

  void add_element(const Type& e)
  {
    elements.push_back(e);
    bounding_box = bounding_box + e.bbox();
  }

  void create_data_structure()
  {
    using boost::bind;
    using boost::make_transform_iterator;

    max_width = CGAL_NTS max BOOST_PREVENT_MACRO_SUBSTITUTION
      (bounding_box.xmax()-bounding_box.xmin(),
       CGAL_NTS max BOOST_PREVENT_MACRO_SUBSTITUTION
       (bounding_box.ymax()-bounding_box.ymin(),
        bounding_box.zmax()-bounding_box.zmin()));
    epsilon = max_width * std::numeric_limits<double>::epsilon();

    std::vector<Interval> intervals;
    for(typename Elements::const_iterator it = elements.begin(),
          end = elements.end(); it != end; ++it)
    {
      intervals.push_back(make_interval(*it));
    }
    tree.make_tree(intervals.begin(),
                   intervals.end());

//     tree.make_tree(make_transform_iterator(elements.begin(),
//                                            bind(&make_interval,_1)),
//                    make_transform_iterator(elements.end(),
//                                            bind(&make_interval,_1)));
    nb_of_elements = elements.size();
    elements.clear();
  }

  const Bbox_3& bbox() const
  {
    return bounding_box;
  }

  double max_length() const
  {
    return max_width;
  }

  typename GT::Iso_cuboid_3 iso_cuboid() const
  {
    return typename GT::Iso_cuboid_3(bbox().xmin(),
                                     bbox().ymin(),
                                     bbox().zmin(),
                                     bbox().xmax(),
                                     bbox().ymax(),
                                     bbox().zmax());
  }

  typename Elements::size_type number_of_elements() const
  {
    return nb_of_elements;
  }

  template <class Type2>
  Object intersection(const Type2& e)
  {
    std::vector<Interval> intervals;

    tree.window_query(std::make_pair(make_pure_interval(e), Type()),
                      std::back_inserter(intervals));
#ifdef CGAL_SURFACE_MESHER_DEBUG_INTERSECTION_DATA_STRUCTURE
    std::cerr << boost::format("intersection percentage: %1$.1f%% query=(%2%, %3%)\n")
      % ( 100. * intervals.size() / number_of_elements() )
      % make_pure_interval(e).first
      % make_pure_interval(e).second;
#endif
    for(typename std::vector<Interval>::const_iterator it = intervals.begin(),
          end = intervals.end(); it != end; ++it)
    {
      Object inter = proper_intersection(it->second, e);
      if(!inter.is_empty())
        return inter;
    }
    return Object();
  }

  template <class Type2>
  std::pair<bool,int> number_of_intersections(const Type2& e)
  {
    int result = 0;

    std::vector<Interval> intervals;

    tree.window_query(std::make_pair(make_pure_interval(e), Type()),
                                     std::back_inserter(intervals));
#ifdef CGAL_SURFACE_MESHER_DEBUG_INTERSECTION_DATA_STRUCTURE
    std::cerr << boost::format("number_of_intersections percentage: %1$.1f%% query=(%2%, %3%)\n")
      % ( 100. * intervals.size() / number_of_elements() )
      % make_pure_interval(e).first
      % make_pure_interval(e).second;
#endif
    for(typename std::vector<Interval>::const_iterator it = intervals.begin(),
          end = intervals.end(); it != end; ++it)
    {
      switch(do_intersect_properly(it->second, e)) {
      case INSIDE:
        ++result;
        break;
      case ON_BOUNDARY:
        return std::make_pair(false, false);
        break;
      case NO_INTERSECTION:
        break;
      }
    }
    return std::make_pair(true,result);
  }
};

} // end namespace CGAL

#endif // CGAL_INTERSECTION_DATA_STRUCTURE_H
