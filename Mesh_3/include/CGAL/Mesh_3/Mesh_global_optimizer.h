// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_GLOBAL_OPTIMIZER_H
#define CGAL_MESH_3_MESH_GLOBAL_OPTIMIZER_H

#ifdef CGAL_MESH_3_VERBOSE
  #define CGAL_MESH_3_OPTIMIZER_VERBOSE 
#endif

#include <CGAL/Timer.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Origin.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/Null_global_optimizer_visitor.h>

#include <vector>
#include <list>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace CGAL {

namespace Mesh_3 {
  
  
template <typename C3T3,
          typename MeshDomain,
          typename MoveFunction,
          typename Visitor_ = Null_global_optimizer_visitor<C3T3> >
class Mesh_global_optimizer
{  
  // Types
  typedef typename C3T3::Triangulation  Tr;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Tr::Point            Point_3;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::Vertex           Vertex;
  
  typedef typename Gt::FT               FT;
  typedef typename Gt::Vector_3         Vector_3;
  
  typedef typename std::vector<Cell_handle>                 Cell_vector;
  typedef typename std::vector<Vertex_handle>               Vertex_vector;
  typedef typename std::set<Vertex_handle>                  Vertex_set;
  typedef std::vector<std::pair<Vertex_handle, Point_3> >   Moves_vector;
  
  typedef typename MoveFunction::Sizing_field Sizing_field;
  
  typedef class C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;
  
  // Visitor class
  // Should define:
  //  - after_compute_moves()
  //  - after_move_points()
  //  - after_rebuild_restricted_delaunay()
  //  - end_of_iteration(int iteration_number)
  typedef Visitor_ Visitor;
    
public:
  /**
   * Constructor
   */
  Mesh_global_optimizer(C3T3& c3t3,
                        const MeshDomain& domain,
                        const FT& freeze_ratio,
                        const FT& convergence_ratio,
                        const MoveFunction move_function = MoveFunction());
  
  /**
   * Launch optimization process
   *
   * @param nb_interations maximum number of iterations
   */
  Mesh_optimization_return_code operator()(int nb_iterations,
                                           Visitor v = Visitor());
  
  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }
  
private:
  /**
   * Returns moves for vertices of set \c moving_vertices
   */
  Moves_vector compute_moves(const Vertex_set& moving_vertices);
  
  /**
   * Returns the move for vertex \c v
   */
  Vector_3 compute_move(const Vertex_handle& v);
  
  /**
   * update big_moves_ vector with new_sq_move value
   */
  void update_big_moves(const FT& new_sq_move);
  
  /**
   * Updates mesh using moves of \c moves vector. Updates moving_vertices with
   * the new set of moving vertices after the move.
   */
  void update_mesh(const Moves_vector& moves,
                   Vertex_set& moving_vertices,
                   Visitor& visitor);
  
  /**
   * Fill sizing field using sizes (avg circumradius) contained in tr_
   */
  void fill_sizing_field();
  
  /**
   * Returns true if convergence is reached
   */
  bool check_convergence() const;
  
  /**
   * Returns the average circumradius length of cells incident to \c v
   */
  FT average_circumradius_length(const Vertex_handle& v) const;
  
  /**
   * Returns the minimum cicumradius length of cells incident to \c v
   */
  FT min_circumradius_sq_length(const Vertex_handle& v) const;
  
  /**
   * Returns the squared circumradius length of cell \c cell
   */
  FT sq_circumradius_length(const Cell_handle& cell,
                            const Vertex_handle& v) const;
  
  /**
   * Returns true if time_limit is reached
   */
  bool is_time_limit_reached() const
  {
    return ( (time_limit() > 0) && (running_time_.time() > time_limit()) );      
  }
  
private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  C3T3& c3t3_;
  Tr& tr_;
  const MeshDomain& domain_;
  FT sq_freeze_ratio_;
  FT convergence_ratio_;
  C3T3_helpers helper_;
  MoveFunction move_function_;
  Sizing_field sizing_field_;
  double time_limit_;
  CGAL::Timer running_time_;
  
  typedef std::list<FT> FT_list;
  FT_list big_moves_;
  
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  mutable FT sum_moves_;
#endif
};
  
  
  
template <typename C3T3, typename Md, typename Mf, typename V_>
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
Mesh_global_optimizer(C3T3& c3t3,
                      const Md& domain,
                      const FT& freeze_ratio,
                      const FT& convergence_ratio,
                      const Mf move_function)
: c3t3_(c3t3)
, tr_(c3t3_.triangulation())
, domain_(domain)
, sq_freeze_ratio_(freeze_ratio*freeze_ratio)
, convergence_ratio_(convergence_ratio)
, helper_(c3t3_,domain_)
, move_function_(move_function)
, sizing_field_(c3t3.triangulation())
, time_limit_(-1)
, running_time_()
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
, sum_moves_(0)
#endif
{
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  std::cerr << "Fill sizing field...";
  CGAL::Timer timer;
  timer.start();
#endif
  
  fill_sizing_field();
  
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  std::cerr << "done (" << timer.time() << "s)\n";
#endif
}


  
template <typename C3T3, typename Md, typename Mf, typename V_>
Mesh_optimization_return_code
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
operator()(int nb_iterations, Visitor visitor)
{
  running_time_.reset();
  running_time_.start();
  
  // Fill set containing moving vertices
  Vertex_set moving_vertices;
  for ( typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
       vit != tr_.finite_vertices_end() ;
       ++vit )
  {
    moving_vertices.insert(vit);
  }

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  double initial_vertices_nb = static_cast<double>(moving_vertices.size());
  double step_begin = running_time_.time();
  
  std::cerr << "Running " << Mf::name() << "-smoothing..." << std::endl;
#endif
  
  // Initialize big moves (stores the largest moves)
  big_moves_.clear();
  std::size_t big_moves_size = (std::max)(std::size_t(1),
                                          moving_vertices.size()/500);
  big_moves_.resize(big_moves_size, FT(0));
  
  // Iterate
  int i = -1;
  while ( ++i < nb_iterations && ! is_time_limit_reached() )
  {
    // Compute move for each vertex
    Moves_vector moves = compute_moves(moving_vertices);
    visitor.after_compute_moves();
    
    // Stop if convergence or time_limit is reached
    if ( check_convergence() || is_time_limit_reached() )
      break;
    
    // Update mesh with those moves
    update_mesh(moves, moving_vertices, visitor);
    visitor.end_of_iteration(i);
    
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
    double moving_vertices_size = static_cast<double>(moving_vertices.size());
    
    std::cerr << boost::format("\r             \r"
                               "end interation %1% (%2$.1f frozen), %3% / %4% (%5%), last step:%6$.2fs, step avg:%7$.2fs, avg large move:%8$.3f          ")
    % i
    % ((1. - moving_vertices_size/initial_vertices_nb)*100.) 
    % moving_vertices_size
    % initial_vertices_nb
    % moves.size() 
    % (running_time_.time() - step_begin)
    % (running_time_.time() / (i+1))
    % sum_moves_;
    
    step_begin = running_time_.time();
#endif
  }
  
  running_time_.stop();
  
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  std::cerr << std::endl;
  
  if ( check_convergence() )
    std::cerr << "Convergence reached" << std::endl;
    
  std::cerr << "Total optimization time: " << running_time_.time()
            << "s" << std::endl << std::endl;
#endif
  
  if ( is_time_limit_reached() )
    return TIME_LIMIT_REACHED;
  
  if ( check_convergence() )
    return CONVERGENCE_REACHED;
  
  return MAX_ITERATION_NUMBER_REACHED;
}

  
template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::Moves_vector
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
compute_moves(const Vertex_set& moving_vertices)
{
  typename Gt::Construct_translated_point_3 translate =
    Gt().construct_translated_point_3_object();
  
  // Store new location of points which have to move
  Moves_vector moves;
  moves.reserve(moving_vertices.size());
  
  // reset worst_move list
  std::fill(big_moves_.begin(),big_moves_.end(),FT(0));
  
  // Get move for each moving vertex
  for ( typename Vertex_set::const_iterator vit = moving_vertices.begin() ;
       vit != moving_vertices.end() ;
       ++vit )
  {
    Vector_3 move = compute_move(*vit);
    if ( CGAL::NULL_VECTOR != move )
    {
      Point_3 new_position = translate((*vit)->point(),move);
      moves.push_back(std::make_pair(*vit,new_position));
    }
    
    // Stop if time_limit_ is reached
    if ( is_time_limit_reached() )
      break;
  }   
  
  return moves;
}
  
  
  
template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::Vector_3
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
compute_move(const Vertex_handle& v)
{    
  typename Gt::Compute_squared_length_3 sq_length =
    Gt().compute_squared_length_3_object();
  
  typename Gt::Construct_vector_3 vector =
    Gt().construct_vector_3_object();
  
  typename Gt::Construct_translated_point_3 translate =
    Gt().construct_translated_point_3_object();
  
  // Get move from move function
  Vector_3 move = move_function_(v, c3t3_, sizing_field_);
  
  // Project surface vertex
  if ( c3t3_.in_dimension(v) == 2 )
  {
    Point_3 new_position = translate(v->point(),move);
    move = vector(v->point(), helper_.project_on_surface(new_position,v));
  }
  
  FT local_sq_size = min_circumradius_sq_length(v);
  if ( FT(0) == local_sq_size )
    return CGAL::NULL_VECTOR;
  
  FT local_move_sq_length = sq_length(move) / local_sq_size;
  
  // Move point only if displacement is big enough w.r.t local size
  if ( local_move_sq_length < sq_freeze_ratio_ )
  {
    return CGAL::NULL_VECTOR;
  }
  
  // Update big moves
  update_big_moves(local_move_sq_length);
  
  return move;
}
  
  
template <typename C3T3, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
update_big_moves(const FT& new_sq_move)
{  
  namespace bl = boost::lambda;
  
  if ( new_sq_move > big_moves_.back() )
  {
    // Remove last value
    big_moves_.pop_back();
    
    // Insert value at the right place
    typename FT_list::iterator pos = 
      std::find_if(big_moves_.begin(), big_moves_.end(), bl::_1 < new_sq_move );
    
    big_moves_.insert(pos, new_sq_move);
  }
}
  
  
template <typename C3T3, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
update_mesh(const Moves_vector& moves,
            Vertex_set& moving_vertices,
            Visitor& visitor)
{ 
  // Cells which have to be updated
  std::set<Cell_handle> outdated_cells;
  
  // Apply moves in triangulation
  for ( typename Moves_vector::const_iterator it = moves.begin() ;
       it != moves.end() ;
       ++it )
  {
    const Vertex_handle& v = it->first;
    const Point_3& new_position = it->second;
    
    // Get size at new position
    if ( Sizing_field::is_vertex_update_needed )
    {
      FT size = sizing_field_(new_position,v);
    
      // Move point
      Vertex_handle new_v = helper_.move_point(v, new_position, outdated_cells);
      
      // Restore size in meshing_info data
      new_v->set_meshing_info(size);
    }
    else
    {
      // Move point
      helper_.move_point(v, new_position, outdated_cells);
    }
    
    // Stop if time_limit_ is reached, here we can't return without rebuilding
    // restricted delaunay
    if ( is_time_limit_reached() )
      break;
  }
  visitor.after_move_points();
  
  // Update c3t3
  moving_vertices.clear();
  helper_.rebuild_restricted_delaunay(outdated_cells.begin(),
                                      outdated_cells.end(),
                                      moving_vertices);
  
  visitor.after_rebuild_restricted_delaunay();
}

  
template <typename C3T3, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
fill_sizing_field()
{
  std::map<Point_3,FT> value_map;
  
  // Fill map with local size
  for(typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin(); 
      vit != tr_.finite_vertices_end();
      ++vit)
  {
    value_map.insert(std::make_pair(vit->point(),
                                    average_circumradius_length(vit)));
  }
  
  // fill sizing field
  sizing_field_.fill(value_map);
}

  
template <typename C3T3, typename Md, typename Mf, typename V_>
bool
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
check_convergence() const
{
  namespace bl = boost::lambda;
  
  FT sum(0);
  for ( typename FT_list::const_iterator
       it = big_moves_.begin(), end = big_moves_.end() ; it != end ; ++it )
  {
    sum += CGAL::sqrt(*it);
  }
  
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  sum_moves_ = sum/big_moves_.size();
#endif
  
  return ( sum/FT(big_moves_.size()) < convergence_ratio_ );
}
  
  
template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::FT
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
average_circumradius_length(const Vertex_handle& v) const
{
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(v, std::back_inserter(incident_cells));
  
  FT sum_len (0);
  unsigned int nb = 0;
  
  for ( typename Cell_vector::iterator cit = incident_cells.begin() ;
       cit != incident_cells.end() ;
       ++cit)
  {
    if ( c3t3_.is_in_complex(*cit) )
    {
      sum_len += CGAL::sqrt(sq_circumradius_length(*cit,v));
      ++nb;
    }
  }
  
  // nb == 0 could happen if there is an isolated point.
  if ( 0 != nb )
  {
    return sum_len/nb;
  }
  else
  {
    // Use outside cells to compute size of point
    for ( typename Cell_vector::iterator cit = incident_cells.begin() ;
         cit != incident_cells.end() ;
         ++cit)
    {
      if ( !tr_.is_infinite(*cit) )
      {
        sum_len += CGAL::sqrt(sq_circumradius_length(*cit,v));
        ++nb;
      }
    }

    CGAL_assertion(nb!=0);
    CGAL_assertion(sum_len!=0);
    return sum_len/nb;
  }
}

  
template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::FT
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
min_circumradius_sq_length(const Vertex_handle& v) const
{
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(v, std::back_inserter(incident_cells));
  
  // Get first cell sq_circumradius_length 
  typename Cell_vector::iterator cit = incident_cells.begin();
  while ( incident_cells.end() != cit && !c3t3_.is_in_complex(*cit) ) { ++cit; }
  
  // if vertex is isolated ...
  if ( incident_cells.end() == cit )
    return FT(0);
  
  // Initialize min
  FT min_sq_len = sq_circumradius_length(*cit++,v);
  
  // Find the minimum value
  for ( ; cit != incident_cells.end() ; ++cit )
  {
    if ( !c3t3_.is_in_complex(*cit) )
      continue;
    
    min_sq_len = (std::min)(min_sq_len,sq_circumradius_length(*cit,v));
  }
  
  return min_sq_len;
}


template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::FT
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
sq_circumradius_length(const Cell_handle& cell, const Vertex_handle& v) const
{
  typename Gt::Compute_squared_distance_3 sq_distance =
    Gt().compute_squared_distance_3_object();
  
  const Point_3 circumcenter = tr_.dual(cell);
  return ( sq_distance(v->point(), circumcenter) );
}
  
} // end namespace Mesh_3

} //namespace CGAL

#endif // CGAL_MESH_3_MESH_GLOBAL_OPTIMIZER_H
