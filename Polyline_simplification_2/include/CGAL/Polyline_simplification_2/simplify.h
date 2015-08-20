// Copyright (c) 2012 Geometry Factory. All rights reserved.
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
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_POLYLINE_SIMPLIFICATION_2_SIMPLIFY_H
#define CGAL_POLYLINE_SIMPLIFICATION_2_SIMPLIFY_H

#include <list>

#include <CGAL/Polyline_simplification_2/Vertex_base_2.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Scaled_squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Hybrid_squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>
#include <CGAL/Polyline_simplification_2/Stop_below_count_threshold.h>
#include <CGAL/Polyline_simplification_2/Stop_above_cost_threshold.h>
#include <CGAL/Modifiable_priority_queue.h>
#include <CGAL/algorithm.h>

// Needed for Polygon_2

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <list>


namespace CGAL {

#ifndef DOXYGEN_RUNNING

template < class CDT >
class Constrained_triangulation_plus_2;


template <class PolygonTraits_2, class Container>
class Polygon_2;

#endif

namespace Polyline_simplification_2 {

template <typename PCT, typename CostFunction, typename StopFunction>
class Polyline_simplification_2
{
public:

  typedef typename PCT::Point Point;
  typedef typename PCT::Constraint_id Constraint_id;
  typedef typename PCT::Constraint_iterator Constraint_iterator;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
  //typedef typename PCT::Points_in_constraint_iterator Points_in_constraint_iterator;
  typedef typename PCT::Vertex_handle Vertex_handle;
  typedef typename PCT::Vertex_circulator Vertex_circulator;

  typedef typename PCT::Geom_traits::FT FT;
  
  PCT& pct;
  CostFunction cost;
  StopFunction stop;
  std::size_t pct_initial_number_of_vertices, number_of_unremovable_vertices;

  
  struct Compare_cost 
  { 
    bool operator() ( Vertices_in_constraint_iterator const& x, 
                      Vertices_in_constraint_iterator const& y ) const 
    { 
      return (*x)->cost() < (*y)->cost(); 
    }
  } ;
  
  struct Id_map : public boost::put_get_helper<std::size_t, Id_map>
  { 
    typedef boost::readable_property_map_tag category;
    typedef std::size_t                      value_type;
    typedef value_type                       reference;
    typedef Vertices_in_constraint_iterator  key_type;
    
    reference operator[] ( key_type const& x ) const { return x.base()->id ; }
  } ;
  
  typedef CGAL::Modifiable_priority_queue<Vertices_in_constraint_iterator,Compare_cost,Id_map> MPQ ;
  
  MPQ* mpq;

  Polyline_simplification_2(PCT& pct, CostFunction cost, StopFunction stop)
    : pct(pct), cost(cost), stop(stop), pct_initial_number_of_vertices(pct.number_of_vertices()), number_of_unremovable_vertices(0)
  {
    int m = initialize_indices();
    initialize_unremovable();
    Compare_cost cc;
    Id_map idm;
    mpq =  new MPQ(m, cc, idm);
    initialize_costs();
  }

  Polyline_simplification_2(PCT& pct, Constraint_id cid, CostFunction cost, StopFunction stop)
    : pct(pct), cost(cost), stop(stop), pct_initial_number_of_vertices(pct.number_of_vertices()), number_of_unremovable_vertices(0)
  {
    int m = initialize_indices(cid);
    initialize_unremovable();
    Compare_cost cc;
    Id_map idm;
    mpq =  new MPQ(m, cc, idm);
    initialize_costs(cid);
  }



  ~Polyline_simplification_2()
  {
    delete mpq;
  }

  void initialize_unremovable()
  {
    std::set<Vertex_handle> vertices;
    Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
    for(; cit!=e; ++cit){
      Constraint_id cid = *cit;
      Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
      (*it)->set_removable(false);
      for(; it != pct.vertices_in_constraint_end(cid); ++it){
        if(vertices.find(*it) != vertices.end()){
          (*it)->set_removable(false);
        } else {
          vertices.insert(*it);
        }
      }
      it = boost::prior(it);
      (*it)->set_removable(false);
    }
  }

  // For all polyline constraints we compute the cost of all unremovable and not removed vertices
  int
  initialize_costs(Constraint_id cid)
  {
    int n=0;
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        ++it){
      if((*it)->is_removable()){
        boost::optional<FT> dist = cost(pct, it);
        if(dist){
          (*it)->set_cost(*dist);
          (*mpq).push(it);
          ++n;
        } else {
          // no need to set the costs as this vertex is not in the priority queue
        } 
      }
    }
    return n;
  }  
  
  void
  initialize_costs()
  {
    int n=0;
    Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
    for(; cit!=e; ++cit){
      n+= initialize_costs(*cit);
    }
  }

  bool
  is_removable(Vertices_in_constraint_iterator it)
  {
    typedef typename PCT::Geom_traits Geom_traits;
    if(! (*it)->is_removable()) {
      return false;
    }
    
    Vertex_handle vh = *it;
    Vertices_in_constraint_iterator u = boost::prior(it);
    Vertex_handle uh = *u;
    Vertices_in_constraint_iterator w = boost::next(it);
    Vertex_handle wh = *w;
    
    typename Geom_traits::Orientation_2 orientation_2 = pct.geom_traits().orientation_2_object();
    CGAL::Orientation o = orientation_2(uh->point(), vh->point(), wh->point());
    if(o == CGAL::COLLINEAR){
      return true;
    }
    if(o == CGAL::LEFT_TURN){
      std::swap(uh,wh);
    }
    
    // uh, vh, wh perform a right turn 
    const Point& up = uh->point();
    const Point& wp = wh->point();
    Vertex_circulator circ = pct.incident_vertices(vh);
    while(circ != uh){
      ++circ;
    }
    ++circ;
    if(circ == wh){
      typename PCT::Edge e;
      CGAL_assertion_code( bool b = ) pct.is_edge(uh,wh,e.first,e.second);
      CGAL_assertion(b);
      return ! pct.is_constrained(e);
    }
    while(circ != wh){
      o = orientation_2(up, circ->point(), wp);
      if(orientation_2(up, wp, circ->point()) != CGAL::RIGHT_TURN){
        return false;
      }
      ++circ;
    }
    return true;
  }


  int
  initialize_indices(Constraint_id cid, int id = 0)
  {
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        ++it){
      it.base()->id = id++;
    }
    return id;
  }
  
  int
  initialize_indices()
  {
    int id = 0;
    Constraint_iterator b = pct.constraints_begin(), e = pct.constraints_end();
    for(; b!=e; ++b){
      id = initialize_indices(*b, id);
    }
    return id;
  }
  
bool
operator()()
{
  if((*mpq).empty()){
      return false;
    }
  Vertices_in_constraint_iterator v = (*mpq).top();
  (*mpq).pop();
  if(stop(pct, *v, (*v)->cost(), pct_initial_number_of_vertices, pct.number_of_vertices())){
    return false;
  }
  if(is_removable(v)){
    Vertices_in_constraint_iterator u = boost::prior(v), w = boost::next(v);
    pct.simplify(v);
    
    if((*u)->is_removable()){
      boost::optional<FT> dist = cost(pct, u);
      if(! dist){
        std::cerr << "undefined cost not handled yet" << std::endl;
      } else {
        (*u)->set_cost(*dist);
        if((*mpq).contains(u)){
        (*mpq).update(u, true);
        }
      }
    }
    
    if((*w)->is_removable()){
      boost::optional<FT> dist = cost(pct, w);
      if(! dist){
        std::cerr << "undefined cost not handled yet" << std::endl;
      } else {
        (*w)->set_cost(*dist);
        if((*mpq).contains(w)){
        (*mpq).update(w, true);
        }

      }
    }
  } else {
    ++number_of_unremovable_vertices;
  }
  return true;
}

  std::size_t
  number_of_removed_vertices() const
  {
    return pct_initial_number_of_vertices - pct.number_of_vertices();
  }

  };

/*!
\ingroup  PkgPolylineSimplification2Functions

Simplifies a single polygon.

\tparam Traits must be a model of `ConstrainedDelaunayTriangulationTraits_2`
\tparam CostFunction must be a model of `PolylineSimplificationCostFunction`.
\tparam StopFunction must be a model of `PolylineSimplificationStopPredicate`

\attention Any \cgal kernel can be used for `Traits`, but as the traits
class is used for internally using a constrained Delaunay triangulation,
it should be a kernel with at least exact predicates.
*/
template <class Traits, class Container, class CostFunction, class StopFunction>
                  CGAL::Polygon_2<Traits,Container>
                  simplify(const CGAL::Polygon_2<Traits,Container>& polygon,
                           CostFunction cost,
                           StopFunction stop)
{
  typedef Traits K;
  typedef typename K::Point_2 Point_2;

  typedef Vertex_base_2< K > Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag> CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT>       PCT;
  typedef typename PCT::Constraint_id Constraint_id;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;

  PCT pct;

  Constraint_id cid = pct.insert_constraint(polygon);

  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(pct, cost, stop);
  while(simplifier()){}

  CGAL::Polygon_2<Traits,Container> result;
  Vertices_in_constraint_iterator beg = pct.vertices_in_constraint_begin(cid);
  Vertices_in_constraint_iterator end = pct.vertices_in_constraint_end(cid);
  for(; beg!=end;){
    Point_2 p = (*beg)->point();
    ++beg;
    if(beg!=end){
      result.push_back(p);
    }
  }
  return result;
}

/*!
\ingroup  PkgPolylineSimplification2Functions

Simplifies an open or closed polyline given as an iterator range of 2D \cgal points.

\tparam PointIterator must be an iterator with value type `CGAL::Kernel::Point_2`.
\tparam CostFunction must be a model of `PolylineSimplificationCostFunction`
\tparam StopFunction must be a model of `PolylineSimplificationStopPredicate`
\tparam PointOutputIterator must be an output iterator to which `CGAL::Kernel::Point_2` can be assigned.
*/
  template <class PointIterator, class CostFunction, class StopFunction, class PointOutputIterator>
  PointOutputIterator
  simplify(PointIterator b, PointIterator e,
           CostFunction cost,
           StopFunction stop,
           PointOutputIterator out,
           bool close = false)
{
  typedef typename std::iterator_traits<PointIterator>::value_type Point_2;
  typedef typename CGAL::Kernel_traits<Point_2>::type K;
  typedef Vertex_base_2< K > Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag> CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT>       PCT;
  typedef typename PCT::Constraint_id Constraint_id;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;

  PCT pct;

  Constraint_id cid = pct.insert_constraint(b,e, close);

  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(pct, cost, stop);
  while(simplifier()){}

  Vertices_in_constraint_iterator beg = pct.vertices_in_constraint_begin(cid);
  Vertices_in_constraint_iterator end = pct.vertices_in_constraint_end(cid);
  for(; beg!=end;){
    Point_2 p = (*beg)->point();
    ++beg;
    if((!close) || (beg!=end)){
      *out++ = p;
    }
  }
  return out;
}


/*!
\ingroup  PkgPolylineSimplification2Functions

Simplifies a single polyline in a triangulation with polylines as constraints. 

\param ct The underlying constrained Delaunay triangulation which embeds the polyline constraints
\param cid The constraint identifier of the polyline constraint to simplify
\param cost The cost function
\param stop The stop function
\param remove_points  If `true` the function \link CGAL::Constrained_triangulation_plus_2::remove_points_without_corresponding_vertex() `ct.remove_points_without_corresponding_vertex()` \endlink is called.
\returns the number of removed vertices
\tparam CDT  must be `CGAL::Constrained_triangulation_plus_2` with a vertex type that
is model of  `PolylineSimplificationVertexBase_2`.
\tparam CostFunction must be a model of `PolylineSimplificationCostFunction`
\tparam StopFunction must be a model of `PolylineSimplificationStopPredicate`
*/
template <class CDT, class CostFunction, class StopFunction>
std::size_t
simplify(CGAL::Constrained_triangulation_plus_2<CDT>& ct,
         typename CGAL::Constrained_triangulation_plus_2<CDT>::Constraint_id cid,
         CostFunction cost,
         StopFunction stop,
         bool remove_points = true)
{
  typedef CGAL::Constrained_triangulation_plus_2<CDT> PCT;
  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(ct, cid, cost, stop);

  while(simplifier()){}
  if(remove_points){
    ct.remove_points_without_corresponding_vertex(cid);
  }
  return simplifier.number_of_removed_vertices();
}

/*!
\ingroup  PkgPolylineSimplification2Functions
Simplifies all polylines in a triangulation with polylines as constraints.
\param ct The underlying constrained Delaunay triangulation which embeds the polyline constraints
\param cost The cost function
\param stop The stop function
\param remove_points If `true` the function \link CGAL::Constrained_triangulation_plus_2::remove_points_without_corresponding_vertex() `ct.remove_points_without_corresponding_vertex()`\endlink is called.
\returns the number of removed vertices
\tparam CDT  must be `CGAL::Constrained_triangulation_plus_2` with a vertex type that
is model of  `PolylineSimplificationVertexBase_2`.
\tparam CostFunction must be a model of `PolylineSimplificationCostFunction`
\tparam StopFunction must be a model of `PolylineSimplificationStopPredicate`
*/

template <class CDT, class CostFunction, class StopFunction>
std::size_t
simplify(CGAL::Constrained_triangulation_plus_2<CDT>& ct,
         CostFunction cost,
         StopFunction stop,
         bool remove_points = true)
{
  typedef CGAL::Constrained_triangulation_plus_2<CDT> PCT;
  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(ct, cost, stop);

  while(simplifier()){}
  if(remove_points){
    ct.remove_points_without_corresponding_vertex();
  }
  return simplifier.number_of_removed_vertices();
}



} // namespace polyline_simplification_2
} // namespace CGAL 
#endif
