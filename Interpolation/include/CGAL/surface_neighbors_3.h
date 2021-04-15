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

#ifndef CGAL_SURFACE_NEIGHBORS_3_H
#define CGAL_SURFACE_NEIGHBORS_3_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Voronoi_intersection_2_traits_3.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Iterator_project.h>

//contains the definition of the Comparator "closer_to_point" and
// the function object Project_vertex_iterator_to_point
#include <CGAL/surface_neighbor_coordinates_3.h>

#include <iterator>
#include <list>
#include <utility>

namespace CGAL {

//without Delaunay filtering
template <class OutputIterator, class InputIterator, class Kernel>
inline
OutputIterator
surface_neighbors_3(InputIterator first, InputIterator beyond,
                    const typename Kernel::Point_3& p,
                    const typename Kernel::Vector_3& normal,
                    OutputIterator out, const Kernel&)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return surface_neighbors_3(first, beyond, p, out, I_gt(p,normal));
}

template <class OutputIterator, class InputIterator, class ITraits>
OutputIterator
surface_neighbors_3(InputIterator first, InputIterator beyond,
                    const typename ITraits::Point_2& p,
                    OutputIterator out, const ITraits& traits)
{
  //definition of the Voronoi intersection triangulation:
  typedef Regular_triangulation_2< ITraits>           I_triangulation;
  typedef typename I_triangulation::Vertex_handle     Vertex_handle;
  typedef typename I_triangulation::Face_handle       Face_handle;
  typedef typename I_triangulation::Locate_type       Locate_type;

  //build Voronoi intersection triangulation:
  I_triangulation it(traits);

  typename ITraits::Construct_weighted_point_2 p2wp =
      it.geom_traits().construct_weighted_point_2_object();
  typename ITraits::Construct_point_2 wp2p =
      it.geom_traits().construct_point_2_object();

  while(first != beyond){
    it.insert(p2wp(*first++));
  }

  const typename ITraits::Weighted_point_2 wp = p2wp(p);

  Locate_type lt;
  int li;
  Face_handle fh = it.locate(wp, lt, li);

  if(lt == I_triangulation::VERTEX){
    *out++ = p;
    return out;
  }

  Vertex_handle vh = it.insert(wp, fh);

  typename I_triangulation::Vertex_circulator vc(it.incident_vertices(vh)),
                                              done(vc);
  do{
    *out++= wp2p(vc->point());
    CGAL_assertion(! it.is_infinite(vc));
  }
  while(vc++ != done);

  return out;
}

//without Delaunay filtering -- certified version:
// a boolean is returned that indicates if a sufficiently large
// neighborhood has been considered so that the
// Voronoi cell of p is not affected by any point outside the smallest
// ball centered on p containing all points in [first,beyond)
template <class OutputIterator, class InputIterator, class Kernel>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first,
                              InputIterator beyond,
                              const typename Kernel::Point_3& p,
                              const typename Kernel::Vector_3& normal,
                              OutputIterator out, const Kernel&)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return surface_neighbors_certified_3(first, beyond, p, out, I_gt(p,normal));
}

//this function takes the radius of the sphere centered on p
// containing the points in [first, beyond] (i.e. the maximal
// distance from p to [first,beyond) as add. parameter:
template <class OutputIterator, class InputIterator, class Kernel>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first,
                              InputIterator beyond,
                              const typename Kernel::Point_3& p,
                              const typename Kernel::Vector_3& normal,
                              const typename Kernel::FT& radius,
                              OutputIterator out,
                              const Kernel& /*K*/)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return surface_neighbors_certified_3(first, beyond, p, radius,
                                       out, I_gt(p,normal));
}

// Versions with instantiated traits class:
template <class OutputIterator, class InputIterator, class ITraits>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first,
                              InputIterator beyond,
                              const typename ITraits::Point_2& p,
                              OutputIterator out,
                              const ITraits& traits)
{
  //find the point in [first,beyond) furthest from p:
  InputIterator furthest = std::max_element(first, beyond,
                                            closer_to_point<ITraits>(p, traits));

  return surface_neighbors_certified_3
      (first, beyond, p,
       traits.compute_squared_distance_2_object()(p,*furthest),
       out, traits);
}

template <class OutputIterator, class InputIterator, class ITraits>
std::pair< OutputIterator, bool >
surface_neighbors_certified_3(InputIterator first,
                              InputIterator beyond,
                              const typename ITraits::Point_2& p,
                              const typename ITraits::FT& radius,
                              OutputIterator out,
                              const ITraits& traits)
{
  //definition of the Voronoi intersection triangulation:
  typedef Regular_triangulation_2< ITraits>      I_triangulation;

  typedef typename I_triangulation::Vertex_handle     Vertex_handle;
  typedef typename I_triangulation::Face_handle       Face_handle;
  typedef typename I_triangulation::Vertex_circulator Vertex_circulator;
  typedef typename I_triangulation::Face_circulator   Face_circulator;
  typedef typename I_triangulation::Locate_type       Locate_type;

  //build Voronoi intersection triangulation:
  I_triangulation it(traits);

  typename ITraits::Construct_weighted_point_2 p2wp =
      it.geom_traits().construct_weighted_point_2_object();
  typename ITraits::Construct_point_2 wp2p =
      it.geom_traits().construct_point_2_object();

  while(first != beyond){
    it.insert(p2wp(*first++));
  }

  const typename ITraits::Weighted_point_2 wp = p2wp(p);

  Locate_type lt;
  int li;
  Face_handle fh = it.locate(wp, lt, li);

  if(lt == I_triangulation::VERTEX){
    *out++ = p;
    return std::make_pair(out,true);
  }
  Vertex_handle vh = it.insert(wp, fh);
  CGAL_assertion(vh->is_valid());

  //determine the furthest distance from p to a vertex of its cell
  bool valid(false);
  Face_circulator fc(it.incident_faces(vh)), fdone(fc);
  do{
    valid = (!it.is_infinite(fc) &&
             (4*radius > traits.compute_squared_distance_2_object()
              (p, it.dual(fc))));
  }while(!valid && ++fc!=fdone);

  //get the neighbor points:
  Vertex_circulator vc(it.incident_vertices(vh)), vdone(vc);
  do{
    *out++ = wp2p(vc->point());
  }while(++vc!=vdone);

  return std::make_pair(out, valid);
}

//using Delaunay triangulation for candidate point filtering:
// => no certification is necessary
template <class Dt, class OutputIterator>
inline
OutputIterator
surface_neighbors_3(const Dt& dt,
                    const typename Dt::Geom_traits::Point_3& p,
                    const typename Dt::Geom_traits::Vector_3& normal,
                    OutputIterator out,
                    typename Dt::Cell_handle start =typename Dt::Cell_handle())
{
  typedef Voronoi_intersection_2_traits_3<typename Dt::Geom_traits> I_gt;
  return surface_neighbors_3(dt, p, out, I_gt(p,normal),start);
}

template <class Dt, class OutputIterator, class ITraits>
OutputIterator
surface_neighbors_3(const Dt& dt,
                    const typename ITraits::Point_2& p,
                    OutputIterator out, const ITraits& traits,
                    typename Dt::Cell_handle start = typename Dt::Cell_handle())
{
  typedef typename ITraits::Point_2       Point_3;

  typedef typename Dt::Cell_handle        Cell_handle;
  typedef typename Dt::Vertex_handle      Vertex_handle;
  typedef typename Dt::Locate_type        Locate_type;

  //the Vertex_handle is, in fact, an iterator over vertex:
  typedef Project_vertex_iterator_to_point< Vertex_handle>    Proj_point;
  typedef Iterator_project<typename std::list< Vertex_handle >::iterator,
                           Proj_point,
                           const Point_3&,
                           const Point_3*,
                           std::ptrdiff_t,
                           std::forward_iterator_tag>         Point_iterator;

  Locate_type lt;
  int li, lj ;
  Cell_handle c = dt.locate(p, lt, li,lj,start);

  //if p is located on a vertex: the only neighbor is found
  if(lt == Dt::VERTEX){
    *out++= (c->vertex(li))->point();
    return out;
  }

  //otherwise get vertices in conflict
  typename std::list< Vertex_handle >  conflict_vertices;
  dt.vertices_on_conflict_zone_boundary(p,c,
                                        std::back_inserter(conflict_vertices));

  for (typename std::list< Vertex_handle >::iterator it = conflict_vertices.begin();
       it != conflict_vertices.end();){
    if(dt.is_infinite(*it)){
      typename std::list< Vertex_handle >::iterator itp = it;
      it++;
      conflict_vertices.erase(itp);
    } else {
      it++;
    }
  }
  return surface_neighbors_3(Point_iterator(conflict_vertices.begin()),
                             Point_iterator(conflict_vertices.end()),
                             p, out, traits);
}

} //namespace CGAL

#endif // CGAL_SURFACE_NEIGHBORS_3_H
