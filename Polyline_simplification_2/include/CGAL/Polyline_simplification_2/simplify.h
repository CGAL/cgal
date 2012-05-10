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

#include <CGAL/Polyline_simplification_2/Scaled_squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Hybrid_squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>
#include <CGAL/Polyline_simplification_2/Stop_below_count_threshold.h>
#include <CGAL/Polyline_simplification_2/Stop_above_cost_threshold.h>
#include <CGAL/Polyline_simplification_2/mark_vertices_unremovable.h>
#include <CGAL/Modifiable_priority_queue.h>
#include <list>

#include <boost/next_prior.hpp>


// Needed for Polygon_2
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyline_constrained_triangulation_2.h>
#include <list>


namespace CGAL {

template < class Tr >
class Polyline_constrained_triangulation_2;


template <class PolygonTraits_2, class Container>
class Polygon_2;

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

  PCT& pct;
  bool keep_points;
  CostFunction cost;
  StopFunction stop;
  std::size_t pct_initial_number_of_vertices, number_of_unremovable_vertices;

  
  struct Compare_cost 
  { 
    bool operator() ( Vertices_in_constraint_iterator const& x, 
                      Vertices_in_constraint_iterator const& y ) const 
    { 
      return (*x)->cost < (*y)->cost; 
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

  Polyline_simplification_2(PCT& pct, CostFunction cost, StopFunction stop, bool keep_points)
    : pct(pct), cost(cost), stop(stop), keep_points(keep_points),  pct_initial_number_of_vertices(pct.number_of_vertices()), number_of_unremovable_vertices(0)
  {
    std::cerr << pct_initial_number_of_vertices << std::endl;
    int m = initialize_indices();
    Compare_cost cc;
    Id_map idm;
    mpq =  new MPQ(m, cc, idm);
    initialize_costs();
  }

  Polyline_simplification_2(PCT& pct, Constraint_id cid, CostFunction cost, StopFunction stop, bool keep_points)
    : pct(pct), cost(cost), stop(stop), keep_points(keep_points),  pct_initial_number_of_vertices(pct.number_of_vertices()), number_of_unremovable_vertices(0)
  {
    int m = initialize_indices(cid);
    Compare_cost cc;
    Id_map idm;
    mpq =  new MPQ(m, cc, idm);
    initialize_costs(cid);
  }



  ~Polyline_simplification_2()
  {
    delete mpq;
  }

  // For all polyline constraints we compute the cost of all non fixed and not removed vertices
  int
  initialize_costs(Constraint_id cid)
  {
    int n=0;
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        ++it){
      if(! (*it)->fixed){
        Vertices_in_constraint_iterator u = boost::prior(it);
        Vertices_in_constraint_iterator w = boost::next(it);
        
        boost::optional<double> dist = cost(pct, u, it, w);
        if(dist){
          (*it)->cost = *dist;
          (*mpq).push(it);
          ++n;
        } else {
          (*it)->cost = (std::numeric_limits<double>::max)();
          std::cerr << "could not compute a cost" << std::endl;
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
  std::cerr << "Initialized cost of " << n << " vertices" << std::endl;
  }

  bool
  is_removable(Vertices_in_constraint_iterator it)
  {
    typedef typename PCT::Geom_traits Geom_traits;
    if((*it)->fixed) {
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
  if(stop(pct, v, (*v)->cost, pct_initial_number_of_vertices, pct.number_of_vertices())){
    return false;
  }
  if(is_removable(v)){
    Vertices_in_constraint_iterator u = boost::prior(v), w = boost::next(v);
    pct.simplify(u,v,w, keep_points);

    if(! (*u)->fixed){
      Vertices_in_constraint_iterator uu = boost::prior(u);
      boost::optional<double> dist = cost(pct, uu,u,w);
      if(! dist){
        std::cerr << "undefined cost not handled yet" << std::endl;
      } else {
        (*u)->cost = *dist;
        if((*mpq).contains(u)){
        (*mpq).update(u, true);
        }
      }
    }
    
    if(! (*w)->fixed){
      Vertices_in_constraint_iterator ww = boost::next(w);
      boost::optional<double> dist = cost(pct, u,w,ww);
      if(! dist){
        std::cerr << "undefined cost not handled yet" << std::endl;
      } else {
        (*w)->cost = *dist;
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


template <class PolygonTraits_2, class Container, class CostFunction, class StopFunction>
                  CGAL::Polygon_2<PolygonTraits_2,Container>
                  simplify(const CGAL::Polygon_2<PolygonTraits_2,Container>& polygon,
                           CostFunction cost,
                           StopFunction stop)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
  typedef CGAL::Exact_predicates_tag                        Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;
  typedef CGAL::Polyline_constrained_triangulation_2<CDT>       PCT;
  typedef PCT::Constraint_id Constraint_id;
  typedef PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;

  PCT pct;

  Constraint_id cid = pct.insert_constraint(polygon);

  mark_vertices_unremovable(pct);
  bool keep_points = false;
  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(pct, cost, stop, keep_points);
  while(simplifier()){}

  CGAL::Polygon_2<PolygonTraits_2,Container> result;
  std::copy(pct.points_in_constraint_begin(cid),
            pct.points_in_constraint_end(cid), std::back_inserter(result));
  return result;
}


template <class Tr, class CostFunction, class StopFunction>
std::size_t
simplify(CGAL::Polyline_constrained_triangulation_2<Tr>& pct,
         typename CGAL::Polyline_constrained_triangulation_2<Tr>::Constraint_id cid,
         CostFunction cost,
         StopFunction stop,
         bool keep_points = false)
{
  typedef CGAL::Polyline_constrained_triangulation_2<Tr> PCT;
  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(pct, cid, cost, stop, keep_points);

  while(simplifier()){}
  if(! keep_points){
    pct.remove_points_from_constraints(cid);
  }
  return simplifier.number_of_removed_vertices();
}


template <class Tr, class CostFunction, class StopFunction>
std::size_t
simplify(CGAL::Polyline_constrained_triangulation_2<Tr>& pct,
         CostFunction cost,
         StopFunction stop,
         bool keep_points = false)
{
  typedef CGAL::Polyline_constrained_triangulation_2<Tr> PCT;
  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(pct, cost, stop, keep_points);

  while(simplifier()){}
  if(! keep_points){
    pct.remove_points_from_constraints();
  }
  std::cerr << "unremovable vertices: " << simplifier.number_of_unremovable_vertices << std::endl;
  std::cerr << "simplify removed " << simplifier.number_of_removed_vertices() << "vertices" << std::endl;
  return simplifier.number_of_removed_vertices();
}

} // namespace polyline_simplification_2
} // namespace CGAL 
#endif
