// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Julia Floetotto
#ifndef CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H
#define CGAL_REGULAR_NEIGHBOR_COORDINATES_2_H

#include <utility>
#include <CGAL/Polygon_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------
//init traits:
template <class Rt, class OutputIterator>
std::pair< OutputIterator, typename Rt::Geom_traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt, 
			       const typename Rt::Geom_traits::
			       Weighted_point& p, 
			       OutputIterator out){
  
  return regular_neighbor_coordinates_2(rt, p, out, 
					rt.geom_traits(),
					typename Rt::Face_handle(NULL));
};
//init start:
template <class Rt, class OutputIterator, class Traits>
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt, 
			       const typename Rt::Geom_traits::
			       Weighted_point& p, 
			       OutputIterator out, const Traits&
			       traits)
{
  
  return regular_neighbor_coordinates_2(rt, p, out,traits, 
					typename Rt::Face_handle(NULL));
};

template <class Rt, class OutputIterator, class Traits, 
  class OutputIteratorVorVertices>
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Traits::Weighted_point& p, 
			       OutputIterator out, OutputIteratorVorVertices
			       vor_vertices, 
			       const Traits& traits){
  return regular_neighbor_coordinates_2(rt, p, out, vor_vertices, traits, 
					typename Rt::Face_handle(NULL));
}
//init vor_vertices to Emptyset_iterator:
template <class Rt, class OutputIterator, class Traits>
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt,
			     const typename Traits::Weighted_point& p, 
			     OutputIterator out, const Traits& traits, 
		 	     typename Rt::Face_handle start){
  
  return regular_neighbor_coordinates_2(rt, p, out, 
					Emptyset_iterator(),
					traits, start);
}

template <class Rt, class OutputIterator, class Traits, 
  class OutputIteratorVorVertices>
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt,
			       const typename Traits::Weighted_point& p, 
			       OutputIterator out, 
			       OutputIteratorVorVertices vor_vertices, 
			       const Traits& traits, 
			       typename Rt::Face_handle start){
  //out: the result of the coordinate computation 
  //vor_vertices: the vertices of the power cell (to avoid recomputation)
  typedef typename Traits::FT            Coord_type;
  typedef typename Traits::Weighted_point  Weighted_point;
  
  typedef typename Rt::Vertex_handle     Vertex_handle;
  typedef typename Rt::Face_handle       Face_handle;
  typedef typename Rt::Edge              Edge;
  typedef typename Rt::Locate_type       Locate_type;
  
  Locate_type lt;
  int li;
  Face_handle fh = rt.locate(p, lt, li, start);

  //the point must lie inside the convex hull:
  CGAL_precondition(lt != Rt::OUTSIDE_AFFINE_HULL && lt !=
		    Rt::OUTSIDE_CONVEX_HULL
		    &&  (!(lt == Rt::EDGE && 
			   (rt.is_infinite(fh) 
			    || rt.is_infinite(fh->neighbor(li))))));
  
  if(lt == Rt::VERTEX){
    //the point must be in conflict:
    CGAL_precondition(rt.power_test(fh->vertex(li)->point(), p) !=
		      ON_NEGATIVE_SIDE); 
    if(rt.power_test(fh->vertex(li)->point(), p) ==ON_ORIENTED_BOUNDARY){
      *out++= std::make_pair(fh->vertex(li)->point(),Coord_type(1));
      return( std::make_pair(out, Coord_type(1)));
    }
  }
  
  std::list<Edge> hole;
  std::list< Vertex_handle > hidden_vertices;
  
  rt.get_boundary_of_conflicts_and_hidden_vertices(p,
						  std::back_inserter(hole),
						  std::back_inserter
						  (hidden_vertices),
						  fh); 
  return 
    regular_neighbor_coordinates_2
    (rt, p, out, vor_vertices, hole.begin(),hole.end(),
     hidden_vertices.begin(), hidden_vertices.end(), traits);
};


template <class Rt, class OutputIterator, class Traits, class
EdgeIterator, class VertexIterator >
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt, 
			       const typename Traits::Weighted_point& p, 
			       OutputIterator out, EdgeIterator
			       hole_begin, EdgeIterator hole_end, 
			       VertexIterator hidden_vertices_begin, 
			       VertexIterator hidden_vertices_end, 
			       const Traits& traits){
  return regular_neighbor_coordinates_2(rt, p,
					out,Emptyset_iterator(), 
					hole_begin, hole_end, 
					hidden_vertices_begin,
					hidden_vertices_end, 
					traits);
}


template <class Rt, class OutputIterator, class Traits, class
EdgeIterator, class VertexIterator , class OutputIteratorVorVertices >
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt, 
			       const typename Traits::Weighted_point& p, 
			       OutputIterator out,  
			       OutputIteratorVorVertices vor_vertices, 
			       EdgeIterator
			       hole_begin, EdgeIterator hole_end, 
			       VertexIterator hidden_vertices_begin, 
			       VertexIterator hidden_vertices_end, 
			       const Traits& traits){
  //precondition: p must lie inside the non-empty hole 
  //               (=^ inside convex hull of neighbors)
  //out: the result of the coordinate computation 
  //vor_vertices: the vertices of the power cell of p (to avoid recomputation)
  CGAL_precondition(rt.dimension()==2);
  
  typedef typename Traits::FT         Coord_type;
  typedef typename Traits::Bare_point      Bare_point;
  typedef typename Traits::Weighted_point  Weighted_point;
  
  typedef typename Rt::Vertex_handle     Vertex_handle;
  typedef typename Rt::Face_circulator   Face_circulator;
  
  //no hole because only (exactly!) one vertex is hidden:
  if(hole_begin==hole_end){
    *out++= std::make_pair((*hidden_vertices_begin)->point(),
			   Coord_type(1));
    ++hidden_vertices_begin;
    CGAL_assertion(hidden_vertices_begin ==hidden_vertices_end); 
    return(std::make_pair(out, Coord_type(1)));
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
      vor[0] = traits.construct_weighted_circumcenter_2_object()
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
	  
	  area += polygon_area_2(vor.begin(), vor.end(), Traits());
	  vor[1] = vor[2];
	}
      //the second Voronoi vertex of the cell of p:
      vor[2] = 
	traits.construct_weighted_circumcenter_2_object()
	(prev->point(),current->point(),p);
      *vor_vertices++= vor[2];
      
      area += polygon_area_2(vor.begin(), vor.end(), Traits());
      *out++= std::make_pair(current->point(),area);
      
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
      area += polygon_area_2(vor.begin(), vor.end(), Traits());
      
      vor[1] = vor[2]; 
      ++fc;
    } 
    
    *out++= std::make_pair((*hidden_vertices_begin)->point(),area);
    area_sum += area;
  }

  return( std::make_pair(out, area_sum));
  
};


/**********************************************************/
//compute the coordinates for a vertex of the triangulation 
// with respect to the other points in the triangulation
template <class Rt, class OutputIterator>
std::pair< OutputIterator, typename Rt::Geom_traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt, 
			       typename Rt::Vertex_handle vh, 
			       OutputIterator out){
  //init the traits class in regular_neighbor_coordinates_2
  // to rt.geom_traits()
  return regular_neighbor_coordinates_2(rt, vh, out, 
					rt.geom_traits());
};
template <class Rt, class OutputIterator, class Traits>
std::pair< OutputIterator, typename Traits::FT > 
regular_neighbor_coordinates_2(const Rt& rt, 
			       typename Rt::Vertex_handle vh, 
			       OutputIterator out, 
			       const Traits& traits){
  //this functions creates a small triangulation of the 
  // incident vertices of this vertex and computes the 
  // natural neighbor coordinates of ch->point() wrt. it.
  typedef typename Rt::Vertex_circulator     Vertex_circulator;
 
  Rt t2;
  Vertex_circulator vc = rt.incident_vertices(vh),
    done(vc);
  do{
    assert(!rt.is_infinite(vc));
    t2.insert(vc->point());
  }
  while(++vc!=done);
    
  return regular_neighbor_coordinates_2(t2, vh->point(), out, 
					traits);
};


//class providing a function object:
template <class Rt, class OutputIterator>
class regular_neighbor_coordinates_2_object 
{
public:
  std::pair< OutputIterator, typename Rt::Geom_traits::FT > 
  operator()(const Rt& rt, 
	     typename Rt::Vertex_handle vh,
	     OutputIterator out){
    return regular_neighbor_coordinates_2(rt, vh, out);
  }
};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif // CGAL_REGULAR_NEIGHBORS_2_H
