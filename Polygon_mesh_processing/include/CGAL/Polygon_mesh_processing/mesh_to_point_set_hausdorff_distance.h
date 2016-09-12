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
//
//
// Author(s)     : Simon Giraudot and Maxime Gimeno

#ifndef MESH_TO_POINT_SET_HAUSDORFF_DISTANCE_H
#define MESH_TO_POINT_SET_HAUSDORFF_DISTANCE_H\

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/convex_hull_incremental_3.h>

template <typename Kernel>
class CPointH
{
private:

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  Point m_point;
  FT m_hausdorff;

public:

  CPointH ()
  {
    m_point = Point (0., 0., 0.);
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

  Point operator() () const { return m_point; }
  Point& operator() () { return m_point; }
  FT hausdorff () const { return m_hausdorff; }

  static Point mid_point (const CPointH& a, const CPointH& b)
  {
    return Point (0.5 * (a().x() + b().x()),
                  0.5 * (a().y() + b().y()),
                  0.5 * (a().z() + b().z()));
  }


};

template <typename Kernel>
class CRefTriangle
{

private:

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  typedef CPointH<Kernel> PointH;
  PointH m_point[3];
  int m_edge;
  FT m_upper_bound;
  FT m_lower_bound;
  Point m_bisector;

public:

  CRefTriangle (const PointH& a, const PointH& b, const PointH& c) // Triangle
  {
    m_point[0] = a;
    m_point[1] = b;
    m_point[2] = c;

    m_edge = -1;
    for (unsigned int i = 0; i < 3; ++ i)
    {
      if (CGAL::angle (m_point[(i+1)%3](), m_point[i](), m_point[(i+2)%3]())
          == CGAL::OBTUSE)
      {
        m_edge = i;
        break;
      }
    }

    if (m_edge == -1)
      m_bisector = CGAL::circumcenter (a(), b(), c());
    else
    {
      m_bisector = PointH::mid_point (m_point[(m_edge+1)%3],
                                      m_point[(m_edge+2)%3]);

      // Point p0 = m_point[(m_edge+1)%3]();
      // Point p1 = m_point[(m_edge+2)%3]();


      // m_bisector = Point ((p0.x () + p1.x ()) / 2.,
      // 		    (p0.y () + p1.y ()) / 2.,
      // 		    (p0.z () + p1.z ()) / 2.);
    }

    m_lower_bound = 0.;
    m_upper_bound = 0.; //std::numeric_limits<FT>::max ();
    for (unsigned int i = 0; i < 3; ++ i)
    {
      if (m_point[i].hausdorff () > m_lower_bound)
        m_lower_bound = m_point[i].hausdorff ();

      FT up = m_point[i].hausdorff ()
        + std::sqrt (CGAL::squared_distance (m_point[i](), m_bisector));

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
  const int edge () const { return m_edge; }

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
};

namespace CGAL{
template <typename Kernel>
class CRefiner
{
private:

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::FT FT;
  typedef CPointH<Kernel> PointH;
  typedef CRefTriangle<Kernel> RefTriangle;
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
    m_upper_bound = std::numeric_limits<FT>::max ();
  }
  ~CRefiner () { }

  bool empty ()
  {
    return m_queue.empty ();
  }
  void reset (FT lower_bound = 0.)
  {
    m_queue = std::priority_queue<RefTriangle> ();
    m_lower_bound = lower_bound;
    m_upper_bound = std::numeric_limits<FT>::max ();
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

    RefTriangle r (PointH (a, std::sqrt (Knn(tree, a, 1).begin()->second)),
                   PointH (b, std::sqrt (Knn(tree, b, 1).begin()->second)),
                   PointH (c, std::sqrt (Knn(tree, c, 1).begin()->second)));

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
    unsigned int before = m_queue.size ();
    //    std::cerr << "Cleaned " << m_queue.size () << " elements to ";
    std::vector<RefTriangle> to_keep;
    while (!(m_queue.empty ()))
    {
      const RefTriangle& current = m_queue.top ();
      if (current.upper_bound () > m_lower_bound)
        to_keep.push_back (current);
      m_queue.pop ();
    }
    //    std::cerr << to_keep.size () << " elements" << std::endl;

    m_queue = std::priority_queue<RefTriangle> ();
    for (auto& r : to_keep)
      m_queue.push (r);

    return (m_queue.size () < before);
  }

  FT refine (FT limit, const Tree& tree, FT upper_bound = 1e30)
  {
    //    std::ofstream f ("rquick.plot");

    unsigned int nb_clean = 0;

    while (uncertainty () > limit && !(m_queue.empty ()))
    {
      if (m_queue.size () > 100000)
      {
        m_queue.top ().print ();
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


      if(CGAL::squared_area(current.points()[0](),
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
      FT hausdorff = std::sqrt (Knn(tree, bisector, 1).begin()->second);
      if (hausdorff > m_lower_bound)
        m_lower_bound = hausdorff;

      if (m_lower_bound > upper_bound)
        return m_upper_bound;

      PointH new_point (bisector, hausdorff);
      int i = current.edge ();
      if (i == -1)
      {
        PointH p0 (current.points()[0]);
        PointH p1 (current.points()[1]);
        PointH p2 (current.points()[2]);

        m_queue.pop ();

        m_queue.push (RefTriangle (new_point, p0, p1));
        m_queue.push (RefTriangle (new_point, p1, p2));
        m_queue.push (RefTriangle (new_point, p2, p0));
      }
      else
      {
        PointH p0 (current.points()[i]);
        PointH p1 (current.points()[(i+1)%3]);
        PointH p2 (current.points()[(i+2)%3]);

        m_queue.pop ();

        m_queue.push (RefTriangle (new_point, p0, p1));
        m_queue.push (RefTriangle (new_point, p0, p2));
      }
    }


    return m_upper_bound;

  }

};
}//CGAL
#endif // MESH_TO_POINT_SET_HAUSDORFF_DISTANCE_H
