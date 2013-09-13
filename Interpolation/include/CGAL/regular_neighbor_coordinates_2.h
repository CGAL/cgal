// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Julia Floetotto

#ifndef CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H
#define CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H

#include <utility>
#include <CGAL/Polygon_2.h>
#include <CGAL/iterator.h>

//for definition of class Project_vertex_output_iterator
#include <CGAL/natural_neighbor_coordinates_2.h>

namespace CGAL {

// in this functions, the traits class is defined via the regular
// triangulation
// see natural_neighbor_coordinates_2 for a proposal for signatures
// that allow to pass the traits class as argument


//the following two functions suppose that
// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
//!!!they are not documented!!!
//init Face_handle start:
// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out)
{
  return regular_neighbor_coordinates_vertex_2(rt, p, out,
					typename Rt::Face_handle());
}

//Face_handle start is known:
// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out,
			       typename Rt::Face_handle start)
{
  return regular_neighbor_coordinates_vertex_2(rt, p, out,
					Emptyset_iterator(), start);
}

//the Voronoi vertices of the power cell are known:
// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator, class OutputIteratorVorVertices>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out,
			       OutputIteratorVorVertices vor_vertices,
			       typename Rt::Face_handle start)
{
  //out: the result of the coordinate computation
  //vor_vertices: the vertices of the power cell (to avoid
  // recomputation)
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

  //the point must lie inside the convex hull
  // sinon return false:
  if(lt == Rt::OUTSIDE_AFFINE_HULL || lt ==
     Rt::OUTSIDE_CONVEX_HULL
     || (lt == Rt::EDGE && (rt.is_infinite(fh)
			    || rt.is_infinite(fh->neighbor(li)))))
    return make_triple(out, Coord_type(1), false);

  if (lt == Rt::VERTEX)
  {
    //the point must be in conflict:
    CGAL_precondition(rt.power_test(fh->vertex(li)->point(), p) !=
		      ON_NEGATIVE_SIDE);
    if (rt.power_test(fh->vertex(li)->point(), p) ==ON_ORIENTED_BOUNDARY)
    {
      *out++= std::make_pair(fh->vertex(li),Coord_type(1));
      return make_triple(out, Coord_type(1), true);
    }
  }

  std::list<Edge> hole;
  std::list< Vertex_handle > hidden_vertices;

  rt.get_boundary_of_conflicts_and_hidden_vertices(p,
						  std::back_inserter(hole),
						  std::back_inserter
						  (hidden_vertices),
						  fh);
  return regular_neighbor_coordinates_vertex_2
            (rt, p, out, vor_vertices, hole.begin(),hole.end(),
             hidden_vertices.begin(), hidden_vertices.end());
}


// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator, class EdgeIterator,
  class VertexIterator >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out, EdgeIterator
			       hole_begin, EdgeIterator hole_end,
			       VertexIterator hidden_vertices_begin,
			       VertexIterator hidden_vertices_end)
{
  return regular_neighbor_coordinates_vertex_2(rt, p,
					out,Emptyset_iterator(),
					hole_begin, hole_end,
					hidden_vertices_begin,
					hidden_vertices_end);
}


// OutputIterator has value type
//        std::pair<Rt::Vertex_handle, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator, class EdgeIterator,
  class VertexIterator , class OutputIteratorVorVertices >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_vertex_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out,
			       OutputIteratorVorVertices vor_vertices,
			       EdgeIterator
			       hole_begin, EdgeIterator hole_end,
			       VertexIterator hidden_vertices_begin,
			       VertexIterator hidden_vertices_end)
{
  //precondition: p must lie inside the non-empty hole
  //               (=^ inside convex hull of neighbors)
  //out: the result of the coordinate computation
  //vor_vertices: the vertices of the power cell of p (to avoid recomputation)
  CGAL_precondition(rt.dimension()==2);

  typedef typename Rt::Geom_traits         Traits;
  typedef typename Traits::FT              Coord_type;
  typedef typename Traits::Bare_point      Bare_point;

  typedef typename Rt::Vertex_handle     Vertex_handle;
  typedef typename Rt::Face_circulator   Face_circulator;

  //no hole because only (exactly!) one vertex is hidden:
  if(hole_begin==hole_end){
    *out++= std::make_pair((*hidden_vertices_begin), Coord_type(1));
    ++hidden_vertices_begin;
    CGAL_assertion(hidden_vertices_begin == hidden_vertices_end);
    return make_triple(out, Coord_type(1), true);
  }

  std::vector<Bare_point>  vor(3);
  Coord_type area_sum(0);

  //determine the last vertex of the hole:
  EdgeIterator hit = hole_end;
  --hit;
  //to start: prev is the "last" vertex of the hole
  // later: prev is the last vertex processed (previously)
  Vertex_handle prev = hit->first->vertex(rt.cw(hit->second));
  hit = hole_begin;
  while(hit != hole_end)
    {
      Coord_type area(0);
      Vertex_handle current = hit->first->vertex(rt.cw(hit->second));

      //a first Voronoi vertex of the cell of p:
      vor[0] = rt.geom_traits().construct_weighted_circumcenter_2_object()
 	(current->point(),
	 hit->first->vertex(rt.ccw(hit->second))->point(), p);
      *vor_vertices++= vor[0];

      //triangulation of the Voronoi subcell:
      //a second vertex as base
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
      vor[2] =
	rt.geom_traits().construct_weighted_circumcenter_2_object()
	(prev->point(),current->point(),p);
      *vor_vertices++= vor[2];

      area += polygon_area_2(vor.begin(), vor.end(), rt.geom_traits());
      *out++= std::make_pair(current,area);

      area_sum += area;

      //update prev and hit:
      prev= current;
      ++hit;
    }

  //get coordinate for hidden vertices
  //                   <=> the area of their Voronoi cell.
  //decomposition of the cell into triangles
  //        vor1: dual of first triangle
  //        vor2, vor 3: duals of two consecutive triangles
  Face_circulator fc, fc_begin;
  for(; hidden_vertices_begin != hidden_vertices_end;
      ++hidden_vertices_begin){
    Coord_type area(0);
    fc_begin = rt.incident_faces(*hidden_vertices_begin);
    vor[0] = rt.dual(fc_begin);
    fc = fc_begin;
    ++fc;
    vor[1] = rt.dual(fc);
    ++fc;
    while(fc != fc_begin){
      vor[2] = rt.dual(fc);
      area += polygon_area_2(vor.begin(), vor.end(), rt.geom_traits());

      vor[1] = vor[2];
      ++fc;
    }

    *out++= std::make_pair((*hidden_vertices_begin),area);
    area_sum += area;
  }

  return make_triple(out, area_sum, true);
}

////////////////////////////////////////////////////////////
//the cast from vertex to point:
// the following functions return an Output_iterator over
// std::pair<Point, FT>
//=> OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
/////////////////////////////////////////////////////////////
//init Face_handle start:
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out)
{
  return regular_neighbor_coordinates_2(rt, p, out,
				        typename Rt::Face_handle());
}

//OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
//Face_handle start is known:
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out,
			       typename Rt::Face_handle start)
{
  return regular_neighbor_coordinates_2(rt, p, out,
				        Emptyset_iterator(), start);
}

//OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
//the Voronoi vertices of the power cell are known:
template <class Rt, class OutputIterator, class OutputIteratorVorVertices>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out,
			       OutputIteratorVorVertices vor_vertices,
			       typename Rt::Face_handle start)
{
  //out: the result of the coordinate computation
  //vor_vertices: the vertices of the power cell (to avoid
  // recomputation)
  Project_vertex_output_iterator<OutputIterator> op(out);

  CGAL_precondition(rt.dimension() == 2);
  
  Triple< Project_vertex_output_iterator<OutputIterator>,
    typename Rt::Geom_traits::FT, bool >  result =
    regular_neighbor_coordinates_vertex_2
    (rt, p, op , vor_vertices, start);
  return make_triple(result.first.base(), result.second, result.third);
}


//OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator, class EdgeIterator,
          class VertexIterator >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out, EdgeIterator
			       hole_begin, EdgeIterator hole_end,
			       VertexIterator hidden_vertices_begin,
			       VertexIterator hidden_vertices_end)
{
   return regular_neighbor_coordinates_2(rt, p,
					out,Emptyset_iterator(),
					hole_begin, hole_end,
					hidden_vertices_begin,
					hidden_vertices_end);
}


//OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator, class EdgeIterator,
  class VertexIterator , class OutputIteratorVorVertices >
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Rt::Weighted_point& p,
			       OutputIterator out,
			       OutputIteratorVorVertices vor_vertices,
			       EdgeIterator hole_begin, EdgeIterator hole_end,
			       VertexIterator hidden_vertices_begin,
			       VertexIterator hidden_vertices_end)
{
  //precondition: p must lie inside the non-empty hole
  //               (=^ inside convex hull of neighbors)
  //out: the result of the coordinate computation
  //vor_vertices: the vertices of the power cell of p
  //(to avoid recomputation)
  Project_vertex_output_iterator<OutputIterator> op(out);

  Triple< Project_vertex_output_iterator<OutputIterator>,
    typename Rt::Geom_traits::FT, bool >  result =
    regular_neighbor_coordinates_vertex_2
    (rt, p, op , vor_vertices, hole_begin,hole_end,
     hidden_vertices_begin, hidden_vertices_end);
  return make_triple(result.first.base(), result.second, result.third);
}

/**********************************************************/
//compute the coordinates for a vertex of the triangulation
// with respect to the other points in the triangulation
//OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator>
Triple< OutputIterator, typename Rt::Geom_traits::FT, bool >
regular_neighbor_coordinates_2(const Rt& rt,
			       typename Rt::Vertex_handle vh,
			       OutputIterator out)
{
  //this functions creates a small triangulation of the
  // incident vertices of this vertex and computes the
  // natural neighbor coordinates of vh->point() wrt. it.
  typedef typename Rt::Vertex_circulator     Vertex_circulator;

  CGAL_precondition(rt.dimension() == 2);
  
  Rt t2;
  Vertex_circulator vc = rt.incident_vertices(vh),
    done(vc);
  do{
    CGAL_assertion(!rt.is_infinite(vc));
    t2.insert(vc->point());
  }
  while(++vc!=done);

  return regular_neighbor_coordinates_2(t2, vh->point(), out);
}


//class providing a function object:
//OutputIterator has value type
//   std::pair< Rt::Geom_traits::Point_2, Rt::Geom_traits::FT>
template <class Rt, class OutputIterator>
class regular_neighbor_coordinates_2_object
{
public:
  Triple< OutputIterator, typename Rt::Geom_traits::FT , bool >
  operator()(const Rt& rt,
	     typename Rt::Vertex_handle vh,
	     OutputIterator out) const
  {
    return regular_neighbor_coordinates_2(rt, vh, out);
  }
};

} //namespace CGAL

#endif // CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H
