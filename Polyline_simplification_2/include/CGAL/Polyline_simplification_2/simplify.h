// Copyright (c) 2012 Geometry Factory. All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_POLYLINE_SIMPLIFICATION_2_SIMPLIFY_H
#define CGAL_POLYLINE_SIMPLIFICATION_2_SIMPLIFY_H

#include <CGAL/license/Polyline_simplification_2.h>

#include <CGAL/disable_warnings.h>

#include <list>
#include <unordered_map>

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

#include <CGAL/Polygon_with_holes_2.h>
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
  typedef typename PCT::Edge Edge;
  typedef typename PCT::Constraint_id Constraint_id;
  typedef typename PCT::Constrained_edges_iterator Constrained_edges_iterator;
  typedef typename PCT::Constraint_iterator Constraint_iterator;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
  typedef typename PCT::Finite_vertices_iterator Finite_vertices_iterator;
  //typedef typename PCT::Points_in_constraint_iterator Points_in_constraint_iterator;
  typedef typename PCT::Vertex_handle Vertex_handle;
    typedef typename PCT::Face_handle Face_handle;
  typedef typename PCT::Vertex_circulator Vertex_circulator;

  typedef typename PCT::Geom_traits::FT FT;

  PCT& pct;
  CostFunction cost;
  StopFunction stop;
  std::size_t pct_initial_number_of_vertices, number_of_unremovable_vertices;

  std::unordered_map<Vertex_handle, std::list<Vertices_in_constraint_iterator> > vertex_to_iterator;

  struct Compare_cost
  {
    bool operator() ( Vertices_in_constraint_iterator const& x,
                      Vertices_in_constraint_iterator const& y ) const
    {
      return (*x)->cost() < (*y)->cost();
    }

    bool operator() (const Vertex_handle& x,const Vertex_handle& y) const
    {
      return x->cost() < y->cost();
    }

  } ;

  struct Id_map
  {
    typedef boost::readable_property_map_tag category;
    typedef std::size_t                      value_type;
    typedef value_type                       reference;
    typedef Vertex_handle                    key_type;


    value_type operator[](const key_type& x) const
    {
      return x->ID;
    }

    friend inline value_type get(const Id_map& m, const key_type k) { return m[k]; }
  } ;

  typedef CGAL::Modifiable_priority_queue<Vertex_handle,Compare_cost,Id_map> MPQ ;

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

  // endpoints of constraints are unremovable
  // vertices which are not endpoint and have != 2 incident constrained edges are unremovable
  void initialize_unremovable()
  {
    Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
    for(; cit!=e; ++cit){
      Constraint_id cid = *cit;
      Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid),
                                      ite = pct.vertices_in_constraint_end(cid);
      (*it)->set_removable(false);
      ++it;
      for(; it != ite; ++it){
        if((boost::next(it) != ite) && (boost::prior(it)== boost::next(it))){
          (*it)->set_removable(false);
        }
      }
      it = boost::prior(it);
      (*it)->set_removable(false);
    }

    std::unordered_map<Vertex_handle, int> degrees;
    for (Constrained_edges_iterator it = pct.constrained_edges_begin(); it != pct.constrained_edges_end(); ++it) {
      Edge e = *it;
      Face_handle fh = e.first;
      int ei = e.second;
      Vertex_handle vh = fh->vertex(pct.cw(ei));
      ++degrees[vh];
      vh = fh->vertex(pct.ccw(ei));
      ++degrees[vh];
    }

    for(Finite_vertices_iterator it = pct.finite_vertices_begin(); it != pct.finite_vertices_end(); ++it){
      if( it->is_removable() && (degrees[it] != 2) ){
        it->set_removable(false);
      }
    }

    cit = pct.constraints_begin(), e = pct.constraints_end();
    for(; cit!=e; ++cit){
      Constraint_id cid = *cit;
      for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
          it != pct.vertices_in_constraint_end(cid);
          ++it){
        if((*it)->is_removable()){
          typename std::unordered_map<Vertex_handle, std::list<Vertices_in_constraint_iterator> >::iterator lit;
          lit = vertex_to_iterator.find(*it);

          if(lit != vertex_to_iterator.end()){
            std::list<Vertices_in_constraint_iterator>& ilist = lit->second;
            if(std::find(ilist.begin(),ilist.end(),it) == ilist.end()){
              ilist.push_back(it);
            }
          }else{
            vertex_to_iterator[*it].push_back(it);
          }
        }
      }
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
          if(! (*mpq).contains(*it)){
              (*mpq).push(*it);
            }
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
        Vertex_handle vh = *it;
        vh->ID = id++;
        //vertex_index_map[vh] = id++;
    }
    return id;
  }

  int
  initialize_indices()
  {
    int id = 0;

    for(Finite_vertices_iterator it =  pct.finite_vertices_begin(); it != pct.finite_vertices_end(); ++it){
      it->ID = id++;
    }
    return id;
  }

bool
operator()()
{
  if((*mpq).empty()){
      return false;
    }
  Vertex_handle v = (*mpq).top();
  (*mpq).pop();
  if(stop(pct, v, v->cost(), pct_initial_number_of_vertices, pct.number_of_vertices())){
    return false;
  }

  Vertices_in_constraint_iterator vit = vertex_to_iterator[v].front();
  if(is_removable(vit)){
    Vertices_in_constraint_iterator u = boost::prior(vit), w = boost::next(vit);
    pct.simplify(vit);

    if((*u)->is_removable()){
      boost::optional<FT> dist = cost(pct, u);
      if(! dist){
        // cost is undefined
        if( mpq->contains(*u) ){
          mpq->erase(*u);
        }
      } else {
        (*u)->set_cost(*dist);
        if(mpq->contains(*u)){
          mpq->update(*u);
        }
        else{
          mpq->push(*u);
        }
      }
    }

    if((*w)->is_removable()){
      boost::optional<FT> dist = cost(pct, w);
      if(! dist){
        // cost is undefined
        if( mpq->contains(*w) ){
          mpq->erase(*w);
        }
      } else {
        (*w)->set_cost(*dist);
        if(mpq->contains(*w)){
          mpq->update(*w);
        }
        else{
          mpq->push(*w);
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
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS,typename internal::Itag<K>::type>  CDT;
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

Simplifies a single polygon with holes.

\tparam Traits must be a model of `ConstrainedDelaunayTriangulationTraits_2`
\tparam CostFunction must be a model of `PolylineSimplificationCostFunction`.
\tparam StopFunction must be a model of `PolylineSimplificationStopPredicate`

\attention Any \cgal kernel can be used for `Traits`, but as the traits
class is used for internally using a constrained Delaunay triangulation,
it should be a kernel with at least exact predicates.
*/
template <class Traits, class Container, class CostFunction, class StopFunction>
CGAL::Polygon_with_holes_2<Traits,Container>
simplify(const CGAL::Polygon_with_holes_2<Traits,Container>& polygon,
         CostFunction cost,
         StopFunction stop)
{
  typedef Traits K;
  typedef typename K::Point_2 Point_2;

  typedef typename CGAL::Polygon_with_holes_2<Traits,Container> Polygon_with_holes_2;
  typedef typename Polygon_with_holes_2::Polygon_2 Polygon_2;

  typedef Vertex_base_2< K > Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS,typename internal::Itag<K>::type>  CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT>       PCT;
  typedef typename PCT::Constraint_id Constraint_id;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;

  PCT pct;

  Constraint_id cid = pct.insert_constraint(polygon.outer_boundary());
  std::vector<Constraint_id> hole_id;
  for(typename Polygon_with_holes_2::Hole_const_iterator it = polygon.holes_begin(); it != polygon.holes_end(); ++it){
     const Polygon_2& hole = *it;
     hole_id.push_back(pct.insert_constraint(hole));
    }

  Polyline_simplification_2<PCT, CostFunction, StopFunction> simplifier(pct, cost, stop);
  while(simplifier()){}

  Polygon_2 result;
  Vertices_in_constraint_iterator beg = pct.vertices_in_constraint_begin(cid);
  Vertices_in_constraint_iterator end = pct.vertices_in_constraint_end(cid);
  for(; beg!=end;){
    Point_2 p = (*beg)->point();
    ++beg;
    if(beg!=end){
      result.push_back(p);
    }
  }
  std::vector<Polygon_2>holes(hole_id.size());
  for(std::size_t i=0; i < hole_id.size(); i++){
    Vertices_in_constraint_iterator beg = pct.vertices_in_constraint_begin(hole_id[i]);
    Vertices_in_constraint_iterator end = pct.vertices_in_constraint_end(hole_id[i]);
    for(; beg!=end;){
      Point_2 p = (*beg)->point();
      ++beg;
      if(beg!=end){
        holes[i].push_back(p);
      }
    }
  }
  return Polygon_with_holes_2(result, holes.begin(), holes.end()) ;
}

/*!
\ingroup  PkgPolylineSimplification2Functions

Simplifies an open or closed polyline given as an iterator range of 2D \cgal points.

\tparam PointIterator must be an iterator with value type `Kernel::Point_2`.
\tparam CostFunction must be a model of `PolylineSimplificationCostFunction`
\tparam StopFunction must be a model of `PolylineSimplificationStopPredicate`
\tparam PointOutputIterator must be an output iterator to which `Kernel::Point_2` can be assigned.
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
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS,typename internal::Itag<K>::type > CDT;
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

#include <CGAL/enable_warnings.h>

#endif
