// Copyright (c) 2016 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Simon Giraudot and Maxime Gimeno

#ifndef CGAL_MESH_TO_POINT_SET_HAUSDORFF_DISTANCE_H
#define CGAL_MESH_TO_POINT_SET_HAUSDORFF_DISTANCE_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>


#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <queue>
#include <iostream>

namespace CGAL{
namespace internal{
template <typename Kernel>
class CPointH
{

public:
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;

  CPointH ()
  {
    m_point = Point (FT(0.), FT(0.), FT(0.));
    m_hausdorff = 0.;
  }
  CPointH (const Point& p, FT h = 0.)
  {
    m_point = p;
    m_hausdorff = h;
  }
  CPointH (const CPointH& p)
  {
    m_point = p ();
    m_hausdorff = p.hausdorff ();
  }

  const Point& operator() () const { return m_point; }
  Point& operator() () { return m_point; }
  FT hausdorff () const { return m_hausdorff; }

  static Point mid_point (const CPointH& a, const CPointH& b)
  {
    return Point (FT(0.5) * (a().x() + b().x()),
                  FT(0.5) * (a().y() + b().y()),
                  FT(0.5) * (a().z() + b().z()));
  }


private:

    Point m_point;
    FT m_hausdorff;
};

template <typename Kernel>
class CRefTriangle
{

private:

  typedef typename Kernel::FT FT;
  typedef CPointH<Kernel> PointH;
  typedef typename PointH::Point Point;
  PointH m_point[3];
  std::size_t m_edge;
  FT m_upper_bound;
  FT m_lower_bound;
  Point m_bisector;

public:

  CRefTriangle (const PointH& a, const PointH& b, const PointH& c) // Triangle
  {
    m_point[0] = a;
    m_point[1] = b;
    m_point[2] = c;

    m_edge = 0;
    
    typename Kernel::Compute_squared_distance_3 squared_distance;
    FT length_max = squared_distance(m_point[1](), m_point[2]());
    FT length1 = squared_distance(m_point[2](), m_point[0]());
    if (length1 > length_max)
      {
        m_edge = 1;
        length_max = length1;
      }
    FT length2 = squared_distance(m_point[0](), m_point[1]());
    if (length2 > length_max)
      m_edge = 2;

    m_bisector = PointH::mid_point (m_point[(m_edge+1)%3],
                                    m_point[(m_edge+2)%3]);

    m_lower_bound = 0.;
    m_upper_bound = 0.;
    for (unsigned int i = 0; i < 3; ++ i)
    {
      if (m_point[i].hausdorff () > m_lower_bound)
        m_lower_bound = m_point[i].hausdorff ();

      FT up = m_point[i].hausdorff ()
        + CGAL::approximate_sqrt (squared_distance (m_point[i](), m_bisector));

      if (up > m_upper_bound)
        m_upper_bound = up;
    }
  }

  CRefTriangle (const CRefTriangle& t)
  {
    m_point[0] = t.points ()[0];
    m_point[1] = t.points ()[1];
    m_point[2] = t.points ()[2];
    m_edge = t.edge ();
    m_lower_bound = t.lower_bound ();
    m_upper_bound = t.upper_bound ();
    m_bisector = t.bisector ();
  }

  FT lower_bound () const
  {
    return m_lower_bound;
  }

  FT upper_bound () const
  {
    return m_upper_bound;
  }

  const Point& bisector () const
  {
    return m_bisector;
  }

  friend bool operator< (const CRefTriangle& a, const CRefTriangle& b)
  {
    return a.upper_bound () < b.upper_bound ();
  }

  const PointH* points () const { return m_point; }
  std::size_t edge () const { return m_edge; }

  #ifdef CGAL_MTPS_HD_DEBUG
  void print () const
  {
    std::cerr << "[Refinement triangle]" << std::endl
              << "   Bounds: " << m_lower_bound << " to "
              << m_upper_bound << std::endl
              << " * " << m_point[0]() << std::endl
              << " * " << m_point[1]() << std::endl
              << " * " << m_point[2]() << std::endl
              << " -> " << m_bisector << std::endl;
  }
  #endif
};
}//internal

template <typename Kernel>
class CRefiner
{
private:

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::FT FT;
  typedef internal::CPointH<Kernel> PointH;
  typedef internal::CRefTriangle<Kernel> RefTriangle;
  typedef CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<Kernel> > Knn;
  typedef typename Knn::Tree Tree;
  std::priority_queue<RefTriangle> m_queue;
  FT m_lower_bound;
  FT m_upper_bound;
  Point m_point;

public:

  CRefiner ()
  {
    m_lower_bound = 0.;
    m_upper_bound = (std::numeric_limits<FT>::max) ();
  }

  bool empty ()
  {
    return m_queue.empty ();
  }
  void reset (FT lower_bound = 0.)
  {
    m_queue = std::priority_queue<RefTriangle> ();
    m_lower_bound = lower_bound;
    m_upper_bound = (std::numeric_limits<FT>::max) ();
  }

  inline FT uncertainty () const
  {
    return m_upper_bound - m_lower_bound;
  }

  inline FT mid () const
  {
    return ((m_upper_bound + m_lower_bound) / 2.);
  }

  inline FT lower_bound () const
  {
    return m_lower_bound;
  }

  inline FT upper_bound () const
  {
    return m_upper_bound;
  }

  inline FT relative_error () const
  {
    return 0.5 * uncertainty () / mid ();
  }


  const Point& hausdorff_point () const { return m_point; }

  bool add (const Point& a, const Point& b, const Point& c, const Tree& tree)
  {

    RefTriangle r (PointH (a, CGAL::approximate_sqrt (Knn(tree, a, 1).begin()->second)),
                   PointH (b, CGAL::approximate_sqrt (Knn(tree, b, 1).begin()->second)),
                   PointH (c, CGAL::approximate_sqrt (Knn(tree, c, 1).begin()->second)));

    if (r.lower_bound () > m_lower_bound)
      {
        m_lower_bound = r.lower_bound ();
      }
    if (r.upper_bound () > m_lower_bound)
      m_queue.push (r);
    return true;
  }
  bool clean_up_queue ()
  {
    std::size_t before = m_queue.size ();
    std::vector<RefTriangle> to_keep;
    while (!(m_queue.empty ()))
    {
      const RefTriangle& current = m_queue.top ();
      if (current.upper_bound () > m_lower_bound)
        to_keep.push_back (current);
      m_queue.pop ();
    }

    m_queue = std::priority_queue<RefTriangle> ();
    BOOST_FOREACH(RefTriangle& r, to_keep)
      m_queue.push (r);

    return (m_queue.size () < before);
  }

  FT refine (FT limit, const Tree& tree, FT upper_bound = 1e30)
  {

    unsigned int nb_clean = 0;

    while (uncertainty () > limit && !(m_queue.empty ()))
    {
      if (m_queue.size () > 100000)
      {
        #ifdef CGAL_MTPS_HD_DEBUG
        m_queue.top ().print ();
        #endif
        ++ nb_clean;
        if (nb_clean > 5)
          return m_upper_bound;
        if (!clean_up_queue ())
          return m_upper_bound;
      }


      const RefTriangle& current = m_queue.top ();

      m_upper_bound = current.upper_bound ();

      if (current.lower_bound () > m_lower_bound)
        m_lower_bound = current.lower_bound ();
      typename Kernel::Compute_squared_area_3 squared_area;

      if(squared_area(current.points()[0](),
                            current.points()[1](),
                            current.points()[2]()
                            ) < 1e-20)
      {
        m_queue.pop();
        continue;
      }
      const Point& bisector = current.bisector ();
      m_point = bisector;
      //squared distance between bisector and its closst point in the mesh
      FT hausdorff = CGAL::approximate_sqrt (Knn(tree, bisector, 1).begin()->second);
      if (hausdorff > m_lower_bound)
        m_lower_bound = hausdorff;

      if (m_lower_bound > upper_bound)
        return m_upper_bound;

      PointH new_point (bisector, hausdorff);
      std::size_t i = current.edge ();
      PointH p0 (current.points()[i]);
      PointH p1 (current.points()[(i+1)%3]);
      PointH p2 (current.points()[(i+2)%3]);

      m_queue.pop ();

      m_queue.push (RefTriangle (new_point, p0, p1));
      m_queue.push (RefTriangle (new_point, p0, p2));
    }


    return m_upper_bound;

  }

};
}//CGAL
#endif // CGAL_MESH_TO_POINT_SET_HAUSDORFF_DISTANCE_H
