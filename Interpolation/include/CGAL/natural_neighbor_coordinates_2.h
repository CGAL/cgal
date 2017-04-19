// Copyright (c) 2003 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Frank Da, Julia Floetotto

#ifndef CGAL_NATURAL_NEIGHBOR_COORDINATES_2_H
#define CGAL_NATURAL_NEIGHBOR_COORDINATES_2_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Iterator_project.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/utility.h>

#include <list>
#include <utility>

namespace CGAL {

// the struct "Project_vertex_output_iterator"
// is used in the (next two) functions
// as well as in regular_neighbor_coordinates_2 and
// in surface_neighbor_coordinates_3
//
//projection of iterator over std::pair <Vertex_handle, T>
//to iterator over std::pair< Point, T>
template < class OutputIterator>
struct Project_vertex_output_iterator
{
  // this class wraps the OutputIterator with value type
  // std::pair<Vertex_handle,T>
  // into an output iterator with value type std::pair<Point, T>
  // Conditions: OutputIterator has value type std::pair<Vertex_handle, T>
  //             and Vertex_handle has a function ->point()
  //             with return type const Point&

  OutputIterator _base;

  //creation:
  Project_vertex_output_iterator(OutputIterator o)
    : _base(o) {}

  OutputIterator base() {return _base;}

  Project_vertex_output_iterator& operator++(){_base++; return *this;}
  Project_vertex_output_iterator& operator++(int){_base++; return *this;}
  Project_vertex_output_iterator& operator*(){return *this;}

  template<class Vertex_pair>
  Project_vertex_output_iterator&
  operator=(const Vertex_pair& x){
    *_base=std::make_pair(x.first->point(), x.second);
    return *this;
  }
};

// The following natural_neighbor_coordinate_2 functions fix the
// traits class to be Dt::Geom_traits. The following signatures could
// be used if one wants to pass a traits class as argument:
// template <class Dt, class OutputIterator, class Traits>
// Triple< OutputIterator, typename Traits::FT, bool >
// natural_neighbor_coordinates_2(const Dt& dt,
//                                const typename Traits::Point_2& p,
//                                OutputIterator out, const Traits& traits,
//                                typename Dt::Face_handle start
//                                = typename Dt::Face_handle())

//template <class Dt, class OutputIterator, class Traits>
//Triple< OutputIterator, typename Traits::FT, bool >
//natural_neighbor_coordinates_2(const Dt& dt,
//                               typename Dt::Vertex_handle vh,
//                               OutputIterator out, const Traits& traits)

//the following two functions suppose that
// OutputIterator has value type
//   std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>
//!!!they are not documented!!!
template < class Dt, class OutputIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_vertex_2(const Dt& dt,
                                      const typename Dt::Geom_traits::Point_2& p,
                                      OutputIterator out, typename Dt::Face_handle start
                                      = typename Dt::Face_handle())
{
  typedef typename Dt::Geom_traits       Traits;
  typedef typename Traits::FT            Coord_type;
  typedef typename Traits::Point_2       Point_2;
  typedef typename Dt::Face_handle       Face_handle;
  typedef typename Dt::Vertex_handle     Vertex_handle;
  typedef typename Dt::Edge              Edge;
  typedef typename Dt::Locate_type       Locate_type;
  typedef typename Traits::Equal_x_2     Equal_x_2;

  CGAL_precondition(dt.dimension() == 2);

  Locate_type lt;
  int li;
  Face_handle fh = dt.locate(p, lt, li, start);

  if (lt == Dt::OUTSIDE_AFFINE_HULL || lt == Dt::OUTSIDE_CONVEX_HULL)
  {
    return make_triple(out, Coord_type(1), false);
  }

  if ((lt == Dt::EDGE &&
       (dt.is_infinite(fh) || dt.is_infinite(fh->neighbor(li)))))
  {
    Vertex_handle v1 = fh->vertex(dt.cw(li));
    Vertex_handle v2 = fh->vertex(dt.ccw(li));

    Point_2 p1(v1->point()), p2(v2->point());

    Coord_type coef1(0);
    Coord_type coef2(0);
    Equal_x_2 equal_x_2;
    if(!equal_x_2(p1,p2))
    {
      coef1 =  (p.x() - p2.x()) / (p1.x() - p2.x());
      coef2 = 1 - coef1;
      *out++ = std::make_pair(v1,coef1);
      *out++ = std::make_pair(v2,coef2);
    } else {
      coef1 = (p.y() - p2.y()) / (p1.y() - p2.y());
      coef2 = 1-coef1;
      *out++ = std::make_pair(v1,coef1);
      *out++ = std::make_pair(v2,coef2);
    }

    return make_triple(out, coef1+coef2, true);
  }

  if (lt == Dt::VERTEX)
  {
    *out++= std::make_pair(fh->vertex(li), Coord_type(1));
    return make_triple(out, Coord_type(1), true);
  }

  std::list<Edge> hole;
  dt.get_boundary_of_conflicts(p, std::back_inserter(hole), fh, false);

  return natural_neighbor_coordinates_vertex_2
      (dt, p, out, hole.begin(), hole.end());
}

//function call if the conflict zone is known:
// OutputIterator has value type
//        std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>
template < class Dt, class OutputIterator, class EdgeIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_vertex_2(const Dt& dt,
                                      const typename Dt::Geom_traits::Point_2& p,
                                      OutputIterator out, EdgeIterator
                                      hole_begin, EdgeIterator hole_end)
{
  CGAL_precondition(dt.dimension()==2);
  //precondition: p must lie inside the hole
  //             (=^ inside convex hull of neighbors)
  typedef typename Dt::Geom_traits       Traits;
  typedef typename Traits::FT            Coord_type;
  typedef typename Traits::Point_2       Point_2;

  typedef typename Dt::Vertex_handle     Vertex_handle;
  typedef typename Dt::Face_circulator   Face_circulator;

  std::vector<Point_2> vor(3);

  Coord_type area_sum(0);
  EdgeIterator hit = hole_end;
  --hit;
  //in the beginning: prev is the "last" vertex of the hole:
  // later: prev is the last vertex processed (previously)
  Vertex_handle prev = hit->first->vertex(dt.cw(hit->second));
  hit = hole_begin;

  while (hit != hole_end)
  {
    Coord_type area(0);
    Vertex_handle current = hit->first->vertex(dt.cw(hit->second));

    vor[0] = dt.geom_traits().construct_circumcenter_2_object()(
               current->point(),
               hit->first->vertex(dt.ccw(hit->second))->point(),
               p);

    Face_circulator fc = dt.incident_faces(current, hit->first);
    ++fc;
    vor[1] = dt.dual(fc);

    while(!fc->has_vertex(prev))
    {
      ++fc;
      vor[2] = dt.dual(fc);
      area += polygon_area_2(vor.begin(), vor.end(), dt.geom_traits());
      vor[1] = vor[2];
    }

    vor[2] = dt.geom_traits().construct_circumcenter_2_object()(prev->point(),
                                                                current->point(),
                                                                p);

    area += polygon_area_2(vor.begin(), vor.end(), dt.geom_traits());
    *out++ = std::make_pair(current,area);
    area_sum += area;

    //update prev and hit:
    prev = current;
    ++hit;
  }
  return make_triple(out, area_sum, true);
}

/////////////////////////////////////////////////////////////
//the cast from vertex to point:
// the following functions return an Output_iterator over
// std::pair<Point, FT>
//=> OutputIterator has value type
//   std::pair< Dt::Geom_traits::Point_2, Dt::Geom_traits::FT>
/////////////////////////////////////////////////////////////
template < class Dt, class OutputIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               OutputIterator out,
                               typename Dt::Face_handle start =
    CGAL_TYPENAME_DEFAULT_ARG Dt::Face_handle() )

{
  CGAL_precondition(dt.dimension() == 2);

  Project_vertex_output_iterator<OutputIterator> op(out);

  Triple<Project_vertex_output_iterator<OutputIterator>,
         typename Dt::Geom_traits::FT, bool > result =
      natural_neighbor_coordinates_vertex_2(dt, p, op, start);

  return make_triple(result.first.base(), result.second, result.third);
}

//OutputIterator has value type
//   std::pair< Dt::Geom_traits::Point_2, Dt::Geom_traits::FT>
//function call if the conflict zone is known:
template < class Dt, class OutputIterator, class EdgeIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               OutputIterator out, EdgeIterator
                               hole_begin, EdgeIterator hole_end)
{
  CGAL_precondition(dt.dimension() == 2);

  Project_vertex_output_iterator<OutputIterator> op(out);

  Triple<Project_vertex_output_iterator<OutputIterator>,
         typename Dt::Geom_traits::FT, bool > result =
      natural_neighbor_coordinates_vertex_2(dt, p, op, hole_begin,hole_end);

  return make_triple(result.first.base(), result.second, result.third);
}

/**********************************************************/
//compute the coordinates for a vertex of the triangulation
// with respect to the other points in the triangulation
//OutputIterator has value type
//   std::pair< Dt::Geom_traits::Point_2, Dt::Geom_traits::FT>
template <class Dt, class OutputIterator>
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               typename Dt::Vertex_handle vh,
                               OutputIterator out)
{
  //this functions creates a small triangulation of the
  // incident vertices of this vertex and computes the
  // natural neighbor coordinates of ch->point() wrt. it.
  typedef typename Dt::Vertex_circulator     Vertex_circulator;

  CGAL_precondition(dt.dimension() == 2);

  Dt t2;
  Vertex_circulator vc = dt.incident_vertices(vh),
      done(vc);
  do{
    CGAL_assertion(!dt.is_infinite(vc));
    t2.insert(vc->point());
  }
  while(++vc!=done);
  return natural_neighbor_coordinates_2(t2, vh->point(), out);
}

//class providing a function object:
//OutputIterator has value type
//   std::pair< Dt::Geom_traits::Point_2, Dt::Geom_traits::FT>
template < class Dt, class OutputIterator >
class natural_neighbor_coordinates_2_object
{
public:
  Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
  operator()(const Dt& dt,
             typename Dt::Vertex_handle vh,
             OutputIterator out) const
  {
    return natural_neighbor_coordinates_2(dt, vh, out);
  }
};

} //namespace CGAL

#endif // CGAL_NATURAL_NEIGHBOR_COORDINATES_2_H
