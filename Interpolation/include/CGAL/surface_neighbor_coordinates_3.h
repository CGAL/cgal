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

// ATTENTION : the surface is supposed to be a closed surface

#ifndef CGAL_SURFACE_NEIGHBOR_COORDINATES_3_H
#define CGAL_SURFACE_NEIGHBOR_COORDINATES_3_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Iterator_project.h>
#include <CGAL/Voronoi_intersection_2_traits_3.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <utility>
#include <vector>

namespace CGAL {

template <class OutputIterator, class InputIterator, class Kernel>
inline
Triple< OutputIterator, typename Kernel::FT, bool >
surface_neighbor_coordinates_3(InputIterator first, InputIterator beyond,
                               const typename Kernel::Point_3& p,
                               const typename Kernel::Vector_3& normal,
                               OutputIterator out,
                               const Kernel&)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return surface_neighbor_coordinates_3(first, beyond, p, out, I_gt(p, normal));
}


template <class OutputIterator, class InputIterator, class ITraits>
Triple< OutputIterator, typename ITraits::FT, bool >
surface_neighbor_coordinates_3(InputIterator first, InputIterator beyond,
                               const typename ITraits::Point_2& p,
                               OutputIterator out,
                               const ITraits& traits)
{
  //definition of the Voronoi intersection triangulation:
  typedef Regular_triangulation_2<ITraits>      I_triangulation;

  //build Voronoi intersection triangulation:
  I_triangulation it(traits);

  typename ITraits::Construct_weighted_point_2 p2wp =
      it.geom_traits().construct_weighted_point_2_object();
  typename ITraits::Construct_point_2 wp2p =
      it.geom_traits().construct_point_2_object();

  while(first != beyond){
    it.insert(p2wp(*first++));
  }

  typedef std::vector<std::pair<typename ITraits::Weighted_point_2,
                                typename ITraits::FT> > WPoint_coordinate_vector;
  WPoint_coordinate_vector wpoint_coords;

  CGAL::Triple<std::back_insert_iterator<WPoint_coordinate_vector>,
               typename ITraits::FT, bool> res =
    regular_neighbor_coordinates_2(it, p2wp(p), std::back_inserter(wpoint_coords));

  // regular_neighbor_coordinates_2 returns weighted_point_2 (in fact, _3)
  // but surface_neighbor_coordinates_3 returns point_3 so the weight must be
  // dropped
  typename WPoint_coordinate_vector::const_iterator wpit = wpoint_coords.begin(),
                                                    wpend = wpoint_coords.end();
  for(; wpit!=wpend; ++wpit)
    *out++ = std::make_pair(wp2p(wpit->first), wpit->second);

  return CGAL::make_triple(out, res.second, res.third);
}

// Without Delaunay filtering but with certification:
// a boolean is returned that indicates if a sufficiently large
// neighborhood has been considered so that the
// Voronoi cell of p is not affected by any point outside the smallest
// ball centered on p containing all points in [first,beyond)
template <class OutputIterator, class InputIterator, class Kernel>
Quadruple< OutputIterator, typename Kernel::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator
                                         first, InputIterator beyond,
                                         const typename Kernel::Point_3& p,
                                         const typename Kernel::Vector_3& normal,
                                         OutputIterator out,
                                         const Kernel& )
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return surface_neighbor_coordinates_certified_3(first, beyond, p, out,
                                                  I_gt(p,normal));
}

//this function takes the radius of the sphere centered on p
// containing the points in [first, beyond] (i.e. the maximal
// distance from p to [first,beyond) as add. parameter:
template <class OutputIterator, class InputIterator, class Kernel>
inline
Quadruple< OutputIterator, typename Kernel::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator first, InputIterator beyond,
                                         const typename Kernel::Point_3& p,
                                         const typename Kernel::Vector_3& normal,
                                         const typename Kernel::FT& radius,
                                         OutputIterator out, const Kernel& )
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return surface_neighbor_coordinates_certified_3(first, beyond, p, radius, out,
                                                  I_gt(p,normal));
}

// FIXME : this should probably be replaced by some kernel functor.
//struct necessary to sort the points by distance to p:
//also used in surface_neighbors_3.h
template <class Traits>
struct closer_to_point
  : public std::less<typename Traits::Point_2>
{
  typedef typename Traits::Point_2   Point_2;

  closer_to_point(const Point_2& _p, const Traits& t)
    : p(_p), traits(t) { }

  bool operator()(const Point_2& q, const Point_2& r) const
  {
    return traits.less_distance_to_point_2_object()(p,q,r);
  }

private:
  Point_2 p;
  Traits traits;
};

// Versions with instantiated traits class:
template <class OutputIterator, class InputIterator, class ITraits>
Quadruple< OutputIterator, typename ITraits::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator first,
                                         InputIterator beyond,
                                         const typename ITraits::Point_2& p,
                                         OutputIterator out,
                                         const ITraits& traits)
{
  //find the point in [first,beyond) furthest from p:
  InputIterator furthest = std::max_element(first, beyond,
                                            closer_to_point<ITraits>(p, traits));

  return surface_neighbor_coordinates_certified_3(first, beyond, p,
           traits.compute_squared_distance_2_object()(p,*furthest),
           out, traits);
}

//with radius(maximal distance from p to [first,beyond)) as
// add. parameter:
template <class OutputIterator, class InputIterator, class ITraits>
Quadruple< OutputIterator, typename ITraits::FT, bool, bool >
surface_neighbor_coordinates_certified_3(InputIterator first,
                                         InputIterator beyond,
                                         const typename ITraits::Point_2& p,
                                         const typename ITraits::FT& radius,
                                         OutputIterator out,
                                         const ITraits& traits)
{
  //definition of the Voronoi intersection triangulation:
  typedef Regular_triangulation_2< ITraits>      I_triangulation;

  //build Voronoi intersection triangulation:
  I_triangulation it(traits);

  typename ITraits::Construct_weighted_point_2 p2wp =
      it.geom_traits().construct_weighted_point_2_object();
  typename ITraits::Construct_point_2 wp2p =
      it.geom_traits().construct_point_2_object();

  while(first != beyond){
    it.insert(p2wp(*first++));
  }

  //collect the Voronoi vertices of the cell of p in order to
  //determine the furthest distance from p to the boundary of its cell
  std::vector< typename ITraits::Point_2 > vor_vertices;

  //unfortunately, there is no function call without Face_handle
  // "start" because this would cause type conflicts because
  // of resembling function signatures (-> default constructor)
  typename ITraits::Weighted_point_2 wp =
    traits.construct_weighted_point_2_object()(p);

  typedef std::vector<std::pair<typename ITraits::Weighted_point_2,
                                typename ITraits::FT> > WPoint_coordinate_vector;
  WPoint_coordinate_vector wpoint_coords;

  Triple< std::back_insert_iterator<WPoint_coordinate_vector>,
          typename ITraits::FT, bool > res =
      regular_neighbor_coordinates_2(it, wp, std::back_inserter(wpoint_coords),
                                     std::back_inserter(vor_vertices),
                                     typename I_triangulation::Face_handle());

  // regular_neighbor_coordinates_2 returns weighted_point_2 (in fact, _3)
  // but surface_neighbor_coordinates_3 returns point_3 so the weight must be
  // dropped
  typename WPoint_coordinate_vector::const_iterator wpit = wpoint_coords.begin(),
                                                    wpend = wpoint_coords.end();
  for(; wpit!=wpend; ++wpit)
    *out++ = std::make_pair(wp2p(wpit->first), wpit->second);

  // if the distance to the furthest sample point is smaller
  // than twice the distance to the furthest vertex, not all neighbors
  // might be found: return false
  typename ITraits::Point_2 furthest =
    *std::max_element(vor_vertices.begin(), vor_vertices.end(),
                      closer_to_point<ITraits>(p,traits));

  if(radius < 4* traits.compute_squared_distance_2_object()(p, furthest))
    return make_quadruple(out, res.second, res.third, false);

  return make_quadruple(out, res.second,res.third, true);
}

// FIXME :
//Sylvain:
//this class should probably be moved to CGAL/function_objects.h
// it is used in the (next two) functions
// it is also used in surface_neighbors_3.h
//
//projection of Vertex_handle (or equiv. Vertices_iterator) to Point
template < class NodeIterator>
struct Project_vertex_iterator_to_point
{
  typedef NodeIterator          argument_type;
  typedef typename std::iterator_traits<NodeIterator>::value_type Node;
  typedef typename Node::Point  Point;
  typedef Point                 result_type;
  Point&       operator()( NodeIterator& x) const { return x->point(); }
  const Point& operator()( const NodeIterator& x) const { return x->point(); }
};

//using Delaunay triangulation for candidate point filtering:
// => no certification is necessary
template <class Dt, class OutputIterator>
inline
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
surface_neighbor_coordinates_3(const Dt& dt,
                               const typename Dt::Geom_traits::Point_3& p,
                               const typename Dt::Geom_traits::Vector_3& normal,
                               OutputIterator out,
                               typename Dt::Cell_handle start = typename Dt::Cell_handle())
{
  typedef Voronoi_intersection_2_traits_3<typename Dt::Geom_traits> I_gt;
  return surface_neighbor_coordinates_3(dt, p, out, I_gt(p,normal), start);
}

template <class Dt, class OutputIterator, class ITraits>
Triple< OutputIterator, typename ITraits::FT, bool >
surface_neighbor_coordinates_3(const Dt& dt,
                               const typename ITraits::Point_2& p,
                               OutputIterator out, const ITraits& traits,
                               typename Dt::Cell_handle start = typename Dt::Cell_handle())
{
  typedef typename ITraits::FT            Coord_type;
  typedef typename ITraits::Point_2       Point_3;

  typedef typename Dt::Cell_handle        Cell_handle;
  typedef typename Dt::Vertex_handle      Vertex_handle;
  typedef typename Dt::Locate_type        Locate_type;

  //the Vertex_handle is, in fact, an iterator over vertex:
  typedef Project_vertex_iterator_to_point< Vertex_handle>   Proj_point;
  typedef Iterator_project<typename std::list< Vertex_handle >::iterator,
                           Proj_point,
                           const Point_3&,
                           const Point_3*,
                           std::ptrdiff_t,
                           std::forward_iterator_tag>  Point_iterator;

  Locate_type lt;
  int li, lj ;
  Cell_handle c = dt.locate(p, lt, li,lj,start);

  //if p is located on a vertex: the only neighbor is found
  if(lt == Dt::VERTEX){
    *out++= std::make_pair(c->vertex(li)->point(), Coord_type(1));
    return make_triple(out, Coord_type(1), true);
  }

  //the candidate points are the points of dt in conflict with p:
  typename std::list< Vertex_handle >  conflict_vertices;
  dt.vertices_on_conflict_zone_boundary(p, c, std::back_inserter(conflict_vertices));

  for (typename std::list< Vertex_handle >::iterator it = conflict_vertices.begin();
       it != conflict_vertices.end();)
  {
    if(dt.is_infinite(*it)){
      typename std::list< Vertex_handle >::iterator itp = it;
      it++;
      conflict_vertices.erase(itp);
    } else {
      it++;
    }
  }
  return surface_neighbor_coordinates_3(Point_iterator(conflict_vertices.begin()),
                                        Point_iterator(conflict_vertices.end()),
                                        p, out, traits);
}

} //namespace CGAL

#endif // CGAL_SURFACE_NEIGHBOR_COORDINATES_3_H
