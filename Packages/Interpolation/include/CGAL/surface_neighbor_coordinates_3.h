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
//
//ATTENTION: 
// the surface is supposed to be a closed surface

#ifndef CGAL_SURFACE_NEIGHBOR_COORDINATES_3_H
#define CGAL_SURFACE_NEIGHBOR_COORDINATES_3_H

#include <utility>

#include <CGAL/Voronoi_intersection_2_traits_3.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>


//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class OutputIterator, class InputIterator, class Kernel>
std::pair< OutputIterator, typename Kernel::FT > 
surface_neighbor_coordinates_3(InputIterator
			       first, InputIterator beyond,
			       const typename Kernel::Point_3& p, 
			       const typename Kernel::Vector_3& normal,
			       OutputIterator out, 
			       const Kernel& K)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return 
    surface_neighbor_coordinates_3(first, beyond, p, out, I_gt(p,normal));
}


template <class OutputIterator, class InputIterator, class ITraits>
std::pair< OutputIterator, typename ITraits::FT > 
surface_neighbor_coordinates_3(InputIterator
			       first, InputIterator beyond,
			       const typename ITraits::Point_2& p,
			       OutputIterator out,  
			       const ITraits& traits)
{
  //definition of the Voronoi intersection triangulation:
  typedef Regular_triangulation_2< ITraits>      I_triangulation;
  
  //build Voronoi intersection triangulation:
  I_triangulation it(traits);
  it.insert(first,beyond);
  
  return 
    regular_neighbor_coordinates_2(it, p, out, traits);
}


//without Delaunay filtering but with certification:
// a boolean is returned that indicates if a sufficiently large 
// neighborhood has been considered so that the
// Voronoi cell of p is not affected by any point outside the smallest 
// ball centered on p containing all points in [first,beyond) 
template <class OutputIterator, class InputIterator, class Kernel>
Triple< OutputIterator, typename Kernel::FT, bool > 
surface_neighbor_coordinates_certified_3(InputIterator
			       first, InputIterator beyond,
			       const typename Kernel::Point_3& p, 
			       const typename Kernel::Vector_3& normal,
			       OutputIterator out, 
			       const Kernel& K)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return 
    surface_neighbor_coordinates_certified_3
    (first, beyond, p, out, I_gt(p,normal));
}
//this function takes the radius of the sphere centered on p 
// containing the points in [first, beyond] (i.e. the maximal 
// distance from p to [first,beyond) as add. parameter: 
template <class OutputIterator, class InputIterator, class Kernel>
Triple< OutputIterator, typename Kernel::FT, bool > 
surface_neighbor_coordinates_certified_3(InputIterator
					 first, InputIterator beyond,
					 const typename Kernel::Point_3& p, 
					 const typename Kernel::Vector_3&
					 normal, 
					 const typename Kernel::FT& 
					 radius,
					 OutputIterator out, 
					 const Kernel& K)
{
  typedef Voronoi_intersection_2_traits_3<Kernel> I_gt;
  return 
    surface_neighbor_coordinates_certified_3
    (first, beyond, p, radius, out, I_gt(p,normal));
}

//struct necessary to sort the points by distance to p:
template <class Traits >
struct closer_to_point : 
  public std::less<typename Traits::Point_2>{
  
  typedef typename Traits::Point_2   Point_2;

  //constructor:
  closer_to_point<Traits>( const Point_2& _p , const Traits& t)
    :p(_p),traits(t){};
  
  bool operator()(const Point_2& q, const Point_2& r)
    {
      return traits.less_distance_to_point_2_object()
	(p,q,r);
    }
private:
  Point_2 p; 
  Traits traits; 
};
// Versions with instantiated traits class:
template <class OutputIterator, class InputIterator, class ITraits>
Triple< OutputIterator, typename ITraits::FT, bool > 
surface_neighbor_coordinates_certified_3(InputIterator
					 first, InputIterator beyond,
					 const typename ITraits::Point_2& p,
					 OutputIterator out,  
					 const ITraits& traits)
{
  //find the point in [first,beyond) furthest from p:
  InputIterator furthest = 
    std::max_element(first, beyond, 
		     closer_to_point<ITraits>(p, traits));
  
  return 
    surface_neighbor_coordinates_certified_3
    (first, beyond, p, 
     traits.compute_squared_distance_2_object()(p,*furthest),
     out, traits);
};

//with radius(maximal distance from p to [first,beyond)) as add. parameter: 
template <class OutputIterator, class InputIterator, class ITraits>
Triple< OutputIterator, typename ITraits::FT, bool > 
surface_neighbor_coordinates_certified_3(InputIterator
					 first, InputIterator beyond,
					 const typename
					 ITraits::Point_2& p,
					 const typename ITraits::FT& radius,
					 OutputIterator out,  
					 const ITraits& traits)
{
  //definition of the Voronoi intersection triangulation:
  typedef Regular_triangulation_2< ITraits>      I_triangulation;
  
  //build Voronoi intersection triangulation:
  I_triangulation it(traits);
  it.insert(first,beyond);
  
  //collect the Voronoi vertices of the cell of p in order to
  //determine the furthest distance from p to the boundary of its cell 
  std::vector<  typename ITraits::Point_2  > vor_vertices;

  std::pair< OutputIterator, typename ITraits::FT > 
    res = regular_neighbor_coordinates_2
    (it, p, out, 
     std::back_inserter(vor_vertices),traits);

   typename ITraits::Point_2 furthest = 
    *std::max_element(vor_vertices.begin(),  vor_vertices.end(), 
		     closer_to_point<ITraits>(p,traits));  
  // if the distance to the furthest sample point is smaller 
  // than twice the distance to the furthest vertex, not all neighbors 
  // might be found: return false
  if(radius < 4* traits.compute_squared_distance_2_object()
     (p, furthest))
    return make_triple(res.first, res.second, false);
  
  return make_triple(res.first, res.second, true);
}

//using Delaunay triangulation for candidate point filtering:
template <class Dt, class OutputIterator>
std::pair< OutputIterator, typename Dt::Geom_traits::FT > 
surface_neighbor_coordinates_3(const Dt& T, 
			       const typename Dt::Geom_traits::Point_3& p,
			       const typename Dt::Geom_traits::Vector_3& 
			       normal, 
			       OutputIterator out, 
			       typename Dt::Cell_handle start 
			       = typename Dt::Cell_handle(NULL)){
  
  typedef Voronoi_intersection_2_traits_3<typename Dt::Geom_traits> I_gt;
  return surface_neighbor_coordinates_3(T, p, out, I_gt(p,normal),start);
  
};

  

template <class Dt, class OutputIterator, class ITraits>
std::pair< OutputIterator, typename ITraits::FT > 
surface_neighbor_coordinates_3(const Dt& T,
			       const typename ITraits::Point_2& p,
			       OutputIterator out, const ITraits& traits,
			       typename Dt::Cell_handle start 
   			       = typename Dt::Cell_handle(NULL)){
  
  typedef typename ITraits::FT            Coord_type;
  typedef typename ITraits::Point_2       Point_3;
  
  
  typedef typename Dt::Cell_handle       Cell_handle;
  typedef typename Dt::Facet             Facet;
  typedef typename Dt::Locate_type       Locate_type;
  
 
 
  Locate_type lt;
  int li, lj ;
  Cell_handle c = T.locate(p, lt, li,lj,start);
  
  if(lt == Dt::VERTEX){
    *out++= std::make_pair(c->vertex(li)->point(),
			   Coord_type(1));
    return( std::make_pair(out, Coord_type(1)));
  }
  
  typename std::vector<Facet> hole;
  T.find_conflicts(p, c, std::back_inserter(hole), 
		   Emptyset_iterator());
  // T.vertices_in_conflict(p, c, std::back_inserter(adjacent));
  
  typename std::set< Point_3, 
    typename Dt::Geom_traits::Less_xyz_3 > adjacent;
  int idx;
  typename std::vector<Facet>::const_iterator hit 
    = hole.begin();
  for(; hit!=hole.end(); ++hit){
    c = hit->first;
    idx = hit->second;
    for(int j=0;j<4;++j)
       if(j!=idx)
	 if(!T.is_infinite(c->vertex(j)))
	   adjacent.insert(c->vertex(j)->point());
  }
  
  return 
    surface_neighbor_coordinates_3
    (adjacent.begin(), adjacent.end(), p, out, traits);
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_NEIGHBORS_3_H
