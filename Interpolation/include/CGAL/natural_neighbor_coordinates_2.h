// Copyright (c) 2003 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Frank Da, Julia Floetotto

#ifndef CGAL_NATURAL_NEIGHBOR_COORDINATES_2_H
#define CGAL_NATURAL_NEIGHBOR_COORDINATES_2_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Interpolation/internal/helpers.h>

#include <CGAL/function_objects.h>
#include <CGAL/is_iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/utility.h>

#include <boost/utility/enable_if.hpp>

#include <iterator>
#include <list>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {

// The following natural_neighbor_coordinate_2 functions fix the
// traits class to be Dt::Geom_traits. The following signatures could
// be used if one wants to pass a traits class as argument:
//
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
// OutputIterator has value type `pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`
//!!!they are not documented!!!
template < class Dt, class OutputIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbors_2(const Dt& dt,
                    const typename Dt::Geom_traits::Point_2& p,
                    OutputIterator out,
                    typename Dt::Face_handle start = typename Dt::Face_handle())
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
    Equal_x_2 equal_x_2;
    if(!equal_x_2(p1,p2))
      coef1 = (p.x() - p2.x()) / (p1.x() - p2.x());
    else
      coef1 = (p.y() - p2.y()) / (p1.y() - p2.y());

    if(coef1 == 0)
    {
      *out++ = std::make_pair(v2, Coord_type(1));
      return make_triple(out, Coord_type(1), true);
    }

    Coord_type coef2 = 1 - coef1;
    if(coef2 == 0)
    {
      *out++ = std::make_pair(v1, Coord_type(1));
      return make_triple(out, Coord_type(1), true);
    }

    *out++ = std::make_pair(v1,coef1);
    *out++ = std::make_pair(v2,coef2);

    return { out, coef1+coef2, true };
  }

  if (lt == Dt::VERTEX)
  {
    *out++= std::make_pair(fh->vertex(li), Coord_type(1));
    return make_triple(out, Coord_type(1), true);
  }

  std::list<Edge> hole;
  dt.get_boundary_of_conflicts(p, std::back_inserter(hole), fh);

  return natural_neighbors_2(dt, p, out, hole.begin(), hole.end());
}


// This function is called when the conflict zone is known.
// OutputIterator has value type `pair<Dt::Vertex_handle, Dt::Geom_traits::FT>`
//
// \pre p must lie inside the hole (=^ inside convex hull of neighbors)
template < class Dt, class OutputIterator, class EdgeIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbors_2(const Dt& dt,
                    const typename Dt::Geom_traits::Point_2& p,
                    OutputIterator out,
                    EdgeIterator hole_begin, EdgeIterator hole_end)
{
  CGAL_precondition(dt.dimension() == 2);

  typedef typename Dt::Geom_traits       Traits;
  typedef typename Traits::FT            Coord_type;
  typedef typename Traits::Point_2       Point_2;

  typedef typename Dt::Vertex_handle     Vertex_handle;
  typedef typename Dt::Face_circulator   Face_circulator;

  std::vector<Point_2> vor(3);

  Coord_type area_sum(0);
  EdgeIterator hit = hole_end;
  --hit;

  // At the beginning, `prev` is the "last" vertex of the hole.
  // Later, `prev` is the last vertex processed (previously).
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

    if(area > 0)
    {
      *out++ = std::make_pair(current,area);
      area_sum += area;
    }

    //update prev and hit:
    prev = current;
    ++hit;
  }
  return make_triple(out, area_sum, true);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// This function is called when the conflict zone is known.
// OutputIterator has value type `pair< Dt::Geom_traits::Point_2, Dt::Geom_traits::FT>`
template < class Dt, class OutputIterator, class OutputFunctor, class EdgeIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               EdgeIterator hole_begin, EdgeIterator hole_end)
{
  CGAL_precondition(dt.dimension() == 2);

  typedef Interpolation::internal::
            Project_vertex_output_iterator<OutputIterator, OutputFunctor> OutputIteratorWithFunctor;
  typedef typename Dt::Geom_traits::FT                                    FT;

  OutputIteratorWithFunctor op(out, fct);
  Triple<OutputIteratorWithFunctor, FT, bool > result = natural_neighbors_2(dt, p, op, hole_begin, hole_end);

  return make_triple(result.first.base(), result.second, result.third);
}


// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template < class Dt, class OutputIterator, class EdgeIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               OutputIterator out,
                               EdgeIterator hole_begin, EdgeIterator hole_end)
{
  typedef typename Dt::Geom_traits::FT                            FT;
  typedef Interpolation::internal::Extract_point_in_pair<Dt, FT>  OutputFunctor;

  return natural_neighbor_coordinates_2(dt, p, out, OutputFunctor(), hole_begin, hole_end);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


// Function with a point. Can take a default starting face.
template < class Dt, class OutputIterator, class OutputFunctor >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               typename Dt::Face_handle start = CGAL_TYPENAME_DEFAULT_ARG Dt::Face_handle(),
                               typename boost::disable_if_c<
                                          is_iterator<OutputFunctor>::value
                                        >::type* = 0)
{
  CGAL_precondition(dt.dimension() == 2);

  typedef Interpolation::internal::
            Project_vertex_output_iterator<OutputIterator, OutputFunctor> OutputIteratorWithFunctor;
  typedef typename Dt::Geom_traits::FT                                    FT;

  OutputIteratorWithFunctor op(out, fct);
  Triple<OutputIteratorWithFunctor, FT, bool > result = natural_neighbors_2(dt, p, op, start);

  return make_triple(result.first.base(), result.second, result.third);
}

// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template < class Dt, class OutputIterator >
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               const typename Dt::Geom_traits::Point_2& p,
                               OutputIterator out,
                               typename Dt::Face_handle start = CGAL_TYPENAME_DEFAULT_ARG Dt::Face_handle())

{
  typedef typename Dt::Geom_traits::FT                            FT;
  typedef Interpolation::internal::Extract_point_in_pair<Dt, FT>  OutputFunctor;

  return natural_neighbor_coordinates_2(dt, p, out, OutputFunctor(), start);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Function with a Vertex_handle
template <class Dt, class OutputIterator, class OutputFunctor>
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               typename Dt::Vertex_handle vh,
                               OutputIterator out,
                               OutputFunctor fct)
{
  // This function creates a small triangulation of the incident vertices of this vertex
  // and computes the natural neighbor coordinates of ch->point() w.r.t. it.
  typedef typename Dt::Vertex_circulator     Vertex_circulator;
  typedef typename Dt::Vertex_handle         Vertex_handle;

  CGAL_precondition(dt.dimension() == 2);

  Dt t2;
  Vertex_circulator vc = dt.incident_vertices(vh), done(vc);

  typedef std::map<Vertex_handle /*t2*/, Vertex_handle/*dt*/> Correspondence_map;
  Correspondence_map cm;

  do
  {
    if(dt.is_infinite(vc))
      continue;

    Vertex_handle t2_vh = t2.insert(vc->point());
    cm[t2_vh] = vc;
  }
  while(++vc != done);

  // Before applying the output functor, we need to switch back from vertices of `t2`
  // to the vertices of `dt`
  typedef typename Dt::Geom_traits::FT                                    FT;
  typedef Interpolation::internal::Pair_mapper<Correspondence_map, FT>    Mapper;
  Mapper pair_mapper(cm);
  CGAL::Unary_compose_1<OutputFunctor, Mapper> composed_fct(fct, pair_mapper);

  return natural_neighbor_coordinates_2(t2, vh->point(), out, composed_fct);
}


// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template <class Dt, class OutputIterator>
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_2(const Dt& dt,
                               typename Dt::Vertex_handle vh,
                               OutputIterator out)
{
  typedef typename Dt::Geom_traits::FT                            FT;
  typedef Interpolation::internal::Extract_point_in_pair<Dt, FT>  OutputFunctor;

  return natural_neighbor_coordinates_2(dt, vh, out, OutputFunctor());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// OutputIterator has value type `std::pair< Dt::Geom_traits::Point_2, Dt::Geom_traits::FT>`
template < class Dt, class OutputIterator, class OutputFunctor >
class natural_neighbor_coordinates_2_object
{
public:
  typedef OutputFunctor                                       Function;

  Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
  operator()(const Dt& dt,
             typename Dt::Vertex_handle vh,
             OutputIterator out,
             OutputFunctor fct) const
  {
    return natural_neighbor_coordinates_2(dt, vh, out, fct);
  }
};

} //namespace CGAL

#endif // CGAL_NATURAL_NEIGHBOR_COORDINATES_2_H
