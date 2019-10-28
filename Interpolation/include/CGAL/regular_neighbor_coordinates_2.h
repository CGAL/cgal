// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H
#define CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Interpolation/internal/helpers.h>

#include <CGAL/is_iterator.h>
#include <CGAL/iterator.h>
#include <CGAL/utility.h>
#include <CGAL/function_objects.h>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <iterator>
#include <list>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {

// In these functions, the traits class is defined via the regular triangulation.
// See natural_neighbor_coordinates_2 for a proposal for signatures
// that allow to pass the traits class as argument.

// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator, class EdgeIterator,
          class VertexIterator , class OutputIteratorVorVertices >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
                                      const typename Rt::Weighted_point& p,
                                      OutputIterator out,
                                      OutputIteratorVorVertices vor_vertices,
                                      EdgeIterator hole_begin,
                                      EdgeIterator hole_end,
                                      VertexIterator hidden_vertices_begin,
                                      VertexIterator hidden_vertices_end)
{
  //precondition: p must lie inside the non-empty hole
  //               (=^ inside convex hull of neighbors)
  //out: the result of the coordinate computation
  //vor_vertices: the vertices of the power cell of p (to avoid recomputation)
  CGAL_precondition(rt.dimension() == 2);

  typedef typename Rt::Geom_traits         Traits;
  typedef typename Traits::FT              Coord_type;
  typedef typename Rt::Bare_point          Bare_point;

  typedef typename Rt::Vertex_handle       Vertex_handle;
  typedef typename Rt::Face_circulator     Face_circulator;

  // no hole
  if(hole_begin == hole_end)
  {
    if(hidden_vertices_begin == hidden_vertices_end)
    {
      // No hole and nothing hidden: the point's weight is too low to appear in the triangulation.
      return make_triple(out, Coord_type(0), true);
    }

    // No hole but some vertices are hidden (there can be only one)
    *out++ = std::make_pair((*hidden_vertices_begin), Coord_type(1));
    ++hidden_vertices_begin;
    CGAL_assertion(hidden_vertices_begin == hidden_vertices_end);
    return make_triple(out, Coord_type(1), true);
  }

  std::vector<Bare_point> vor(3);
  Coord_type area_sum(0);

  //determine the last vertex of the hole:
  EdgeIterator hit = hole_end;
  --hit;

  // to start: prev is the "last" vertex of the hole
  // later: prev is the last vertex processed (previously)
  Vertex_handle prev = hit->first->vertex(rt.cw(hit->second));
  hit = hole_begin;
  while(hit != hole_end)
  {
    Coord_type area(0);
    Vertex_handle current = hit->first->vertex(rt.cw(hit->second));

    // a first Voronoi vertex of the cell of p:
    vor[0] = rt.geom_traits().construct_weighted_circumcenter_2_object()(
               current->point(), hit->first->vertex(rt.ccw(hit->second))->point(), p);
    *vor_vertices++= vor[0];

    // triangulation of the Voronoi subcell:
    // a second vertex as base
    Face_circulator fc = rt.incident_faces(current, hit->first);
    ++fc;
    vor[1] = rt.dual(fc);

    // iteration over all other "old" Voronoi vertices
    while(!fc->has_vertex(prev))
    {
      ++fc;
      vor[2] = rt.dual(fc);

      area += polygon_area_2(vor.begin(), vor.end(), rt.geom_traits());
      vor[1] = vor[2];
    }

    //the second Voronoi vertex of the cell of p:
    vor[2] = rt.geom_traits().construct_weighted_circumcenter_2_object()(
               prev->point(),current->point(), p);
    *vor_vertices++ = vor[2];

    area += polygon_area_2(vor.begin(), vor.end(), rt.geom_traits());

    if(area > 0)
    {
      *out++= std::make_pair(current, area);
      area_sum += area;
    }

    //update prev and hit:
    prev = current;
    ++hit;
  }

  // get coordinate for hidden vertices <=> the area of their Voronoi cell.
  // decomposition of the cell into triangles
  //         vor1: dual of first triangle
  //         vor2, vor 3: duals of two consecutive triangles
  Face_circulator fc, fc_begin;
  for(; hidden_vertices_begin != hidden_vertices_end; ++hidden_vertices_begin)
  {
    Coord_type area(0);
    fc_begin = rt.incident_faces(*hidden_vertices_begin);
    vor[0] = rt.dual(fc_begin);
    fc = fc_begin;
    ++fc;
    vor[1] = rt.dual(fc);
    ++fc;
    while(fc != fc_begin)
    {
      vor[2] = rt.dual(fc);
      area += polygon_area_2(vor.begin(), vor.end(), rt.geom_traits());

      vor[1] = vor[2];
      ++fc;
    }

    if(area > 0)
    {
      *out++ = std::make_pair((*hidden_vertices_begin), area);
      area_sum += area;
    }
  }

  return make_triple(out, area_sum, true);
}


// OutputIterator has value type `pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator, class EdgeIterator,
          class VertexIterator >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
                                      const typename Rt::Weighted_point& p,
                                      OutputIterator out,
                                      EdgeIterator hole_begin,
                                      EdgeIterator hole_end,
                                      VertexIterator hidden_vertices_begin,
                                      VertexIterator hidden_vertices_end)
{
  return regular_neighbor_coordinates_vertex_2(rt, p, out, Emptyset_iterator(),
                                               hole_begin, hole_end,
                                               hidden_vertices_begin,
                                               hidden_vertices_end);
}


// vor_vertices: the vertices of the power cell (to avoid recomputation)
// OutputIterator has value type `pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator, class OutputIteratorVorVertices>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
                                      const typename Rt::Weighted_point& p,
                                      OutputIterator out,
                                      OutputIteratorVorVertices vor_vertices,
                                      typename Rt::Face_handle start)
{
  typedef typename Rt::Geom_traits        Traits;
  typedef typename Traits::FT             Coord_type;

  typedef typename Rt::Vertex_handle      Vertex_handle;
  typedef typename Rt::Face_handle        Face_handle;
  typedef typename Rt::Edge               Edge;
  typedef typename Rt::Locate_type        Locate_type;

  CGAL_precondition(rt.dimension() == 2);

  Locate_type lt;
  int li;
  Face_handle fh = rt.locate(p, lt, li, start);

  // the point must lie inside the convex hull, otherwise return false:
  if(lt == Rt::OUTSIDE_AFFINE_HULL || lt == Rt::OUTSIDE_CONVEX_HULL
     || (lt == Rt::EDGE && (rt.is_infinite(fh) || rt.is_infinite(fh->neighbor(li)))))
    return make_triple(out, Coord_type(1), false);

  if(lt == Rt::VERTEX)
  {
    //the point must be in conflict:
    CGAL_precondition(rt.power_test(fh->vertex(li)->point(), p) != ON_NEGATIVE_SIDE);
    if (rt.power_test(fh->vertex(li)->point(), p) == ON_ORIENTED_BOUNDARY)
    {
      *out++ = std::make_pair(fh->vertex(li), Coord_type(1));
      return make_triple(out, Coord_type(1), true);
    }
  }

  std::list<Edge> hole;
  std::list<Vertex_handle> hidden_vertices;
  rt.get_boundary_of_conflicts_and_hidden_vertices(p,
                                                   std::back_inserter(hole),
                                                   std::back_inserter(hidden_vertices),
                                                   fh);

  return regular_neighbor_coordinates_vertex_2(rt, p, out, vor_vertices,
                                               hole.begin(),hole.end(),
                                               hidden_vertices.begin(),
                                               hidden_vertices.end());
}


// OutputIterator has value type `pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
                                      const typename Rt::Weighted_point& p,
                                      OutputIterator out,
                                      typename Rt::Face_handle start)
{
  return regular_neighbor_coordinates_vertex_2(rt, p, out, Emptyset_iterator(), start);
}


// OutputIterator has value type `pair<Rt::Vertex_handle, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
                                      const typename Rt::Weighted_point& p,
                                      OutputIterator out)
{
  return regular_neighbor_coordinates_vertex_2(rt, p, out, typename Rt::Face_handle());
}

////////////////////////////////////////////////////////////
//the cast from vertex to point:
// the following functions return an Output_iterator over
// std::pair<Point, FT>
//=> OutputIterator has value type
//   std::pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>
/////////////////////////////////////////////////////////////

// OutputIterator has value type `pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>`
// out: the result of the coordinate computation
// vor_vertices: the vertices of the power cell (to avoid recomputation)
template < class Rt, class OutputIterator, class OutputFunctor, class OutputIteratorVorVertices >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               OutputIteratorVorVertices vor_vertices,
                               typename Rt::Face_handle start)
{
  typedef Interpolation::internal::
            Project_vertex_output_iterator<OutputIterator, OutputFunctor> OutputIteratorWithFunctor;
  typedef typename Rt::Geom_traits::FT                                    FT;

  OutputIteratorWithFunctor op(out, fct);

  CGAL_precondition(rt.dimension() == 2);

  Triple< OutputIteratorWithFunctor, FT, bool > result =
    regular_neighbor_coordinates_vertex_2(rt, p, op, vor_vertices, start);

  return make_triple(result.first.base(), result.second, result.third);
}


// OutputIteratorVorVertices and start are given but not OutputFunctor.
template <class Rt, class OutputIterator, class OutputIteratorVorVertices>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputIteratorVorVertices vor_vertices,
                               typename Rt::Face_handle start,
                               typename boost::enable_if_c<
                                          is_iterator<OutputIteratorVorVertices>::value
                                        >::type* = 0)
{
  // Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
  typedef typename Rt::Geom_traits::FT                            FT;
  typedef Interpolation::internal::Extract_point_in_pair<Rt, FT>  OutputFunctor;

  return regular_neighbor_coordinates_2(rt, p, out, OutputFunctor(), vor_vertices, start);
}


// OutputFunctor and start are given but not OutputIteratorVorVertices
template < class Rt, class OutputIterator, class OutputFunctor >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               typename Rt::Face_handle start,
                               typename boost::disable_if_c<
                                          is_iterator<OutputFunctor>::value
                                        >::type* = 0)
{
  return regular_neighbor_coordinates_2(rt, p, out, fct, Emptyset_iterator(), start);
}


// Only the OutputFunctor is given.
template < class Rt, class OutputIterator, class OutputFunctor >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               typename boost::disable_if_c<
                                 boost::is_convertible<OutputFunctor,
                                                       typename Rt::Face_handle>::value
                               >::type* = 0)
{
  return regular_neighbor_coordinates_2(rt, p, out, fct, typename Rt::Face_handle());
}


// Only the starting face is given.
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               typename Rt::Face_handle start)
{
  return regular_neighbor_coordinates_2(rt, p, out, Emptyset_iterator(), start);
}


// Nothing is given.
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out)
{
  return regular_neighbor_coordinates_2(rt, p, out, typename Rt::Face_handle());
}


////////////////////////////////////////////////////////////
// Functions below use the conflict region
////////////////////////////////////////////////////////////

// OutputIterator has value type `pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>`
// vor_vertices are the vertices of the power cell of p (to avoid recomputation)
//
// \pre p must lie inside the non-empty hole (=^ inside convex hull of neighbors)
template < class Rt, class OutputIterator, class OutputFunctor,
           class EdgeIterator, class VertexIterator, class OutputIteratorVorVertices >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               OutputIteratorVorVertices vor_vertices,
                               EdgeIterator hole_begin,
                               EdgeIterator hole_end,
                               VertexIterator hidden_vertices_begin,
                               VertexIterator hidden_vertices_end)
{
  typedef Interpolation::internal::
            Project_vertex_output_iterator<OutputIterator, OutputFunctor> OutputIteratorWithFunctor;
  typedef typename Rt::Geom_traits::FT                                    FT;

  OutputIteratorWithFunctor op(out, fct);

  Triple<OutputIteratorWithFunctor, FT, bool> result =
    regular_neighbor_coordinates_vertex_2(rt, p, op, vor_vertices,
                                          hole_begin, hole_end,
                                          hidden_vertices_begin,
                                          hidden_vertices_end);

  return make_triple(result.first.base(), result.second, result.third);
}


// Same as above but without OutputIteratorVorVertices.
// OutputIterator has value type `pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator, class OutputFunctor,
          class EdgeIterator, class VertexIterator >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputFunctor fct,
                               EdgeIterator hole_begin,
                               EdgeIterator hole_end,
                               VertexIterator hidden_vertices_begin,
                               VertexIterator hidden_vertices_end,
                               typename boost::disable_if_c<
                                          is_iterator<OutputFunctor>::value
                                        >::type* = 0)
{
   return regular_neighbor_coordinates_2(rt, p, out, fct, Emptyset_iterator(),
                                         hole_begin, hole_end,
                                         hidden_vertices_begin,
                                         hidden_vertices_end);
}


// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template < class Rt, class OutputIterator, class OutputIteratorVorVertices,
           class EdgeIterator, class VertexIterator >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               OutputIteratorVorVertices vor_vertices,
                               EdgeIterator hole_begin,
                               EdgeIterator hole_end,
                               VertexIterator hidden_vertices_begin,
                               VertexIterator hidden_vertices_end,
                               typename boost::enable_if_c<
                                          is_iterator<OutputIteratorVorVertices>::value
                                        >::type* = 0)
{
  typedef typename Rt::Geom_traits::FT                            FT;
  typedef Interpolation::internal::Extract_point_in_pair<Rt, FT>  OutputFunctor;

  return regular_neighbor_coordinates_2(rt, p, out, OutputFunctor(),
                                        vor_vertices, hole_begin, hole_end,
                                        hidden_vertices_begin, hidden_vertices_end);
}


// Same as above but without OutputFunctor nor OutputIteratorVorVertices
// OutputIterator has value type `pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator, class EdgeIterator, class VertexIterator >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               const typename Rt::Weighted_point& p,
                               OutputIterator out,
                               EdgeIterator hole_begin,
                               EdgeIterator hole_end,
                               VertexIterator hidden_vertices_begin,
                               VertexIterator hidden_vertices_end)
{
   return regular_neighbor_coordinates_2(rt, p, out, Emptyset_iterator(),
                                         hole_begin, hole_end,
                                         hidden_vertices_begin,
                                         hidden_vertices_end);
}


/**********************************************************/
// Compute the coordinates for a vertex of the triangulation
// with respect to the other points in the triangulation
// OutputIterator has value type `pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>`
template <class Rt, class OutputIterator, class OutputFunctor>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               typename Rt::Vertex_handle vh,
                               OutputIterator out,
                               OutputFunctor fct)
{
  // This function creates a small triangulation of the incident vertices of this vertex
  // and computes the natural neighbor coordinates of ch->point() w.r.t. it.
  typedef typename Rt::Vertex_circulator     Vertex_circulator;
  typedef typename Rt::Vertex_handle         Vertex_handle;

  CGAL_precondition(rt.dimension() == 2);

  Rt t2;
  Vertex_circulator vc = rt.incident_vertices(vh), done(vc);

  typedef std::map<Vertex_handle /*t2*/, Vertex_handle/*dt*/> Correspondence_map;
  Correspondence_map cm;

  do
  {
    if(rt.is_infinite(vc))
      continue;

    Vertex_handle t2_vh = t2.insert(vc->point());
    cm[t2_vh] = vc;
  }
  while(++vc != done);

  // Before applying the output functor, we need to switch back from vertices of `t2`
  // to the vertices of `rt`
  typedef typename Rt::Geom_traits::FT                                    FT;
  typedef Interpolation::internal::Pair_mapper<Correspondence_map, FT>    Mapper;
  Mapper pair_mapper(cm);
  CGAL::Unary_compose_1<OutputFunctor, Mapper> composed_fct(fct, pair_mapper);

  return regular_neighbor_coordinates_2(t2, vh->point(), out, composed_fct);
}


// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
                               typename Rt::Vertex_handle vh,
                               OutputIterator out)
{
  typedef typename Rt::Geom_traits::FT                            FT;
  typedef Interpolation::internal::Extract_point_in_pair<Rt, FT>  OutputFunctor;

  return regular_neighbor_coordinates_2(rt, vh, out, OutputFunctor());
}


// OutputIterator has value type `pair< Rt::Geom_traits::Weighted_point_2, Rt::Geom_traits::FT>`
template < class Rt, class OutputIterator, class OutputFunctor >
class regular_neighbor_coordinates_2_object
{
public:
  typedef OutputFunctor                                       Function;

  Triple< OutputIterator, typename Rt::Geom_traits::FT , bool >
  operator()(const Rt& rt,
             typename Rt::Vertex_handle vh,
             OutputIterator out,
             OutputFunctor fct) const
  {
    return regular_neighbor_coordinates_2(rt, vh, out, fct);
  }
};

} //namespace CGAL

#endif // CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H
