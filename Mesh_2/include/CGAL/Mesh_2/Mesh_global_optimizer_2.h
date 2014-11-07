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
// Author(s)     : Jane Tournois, Raul Gallegos, Stéphane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_2_H
#define CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_2_H

#ifdef CGAL_MESH_2_VERBOSE
  #define CGAL_MESH_2_OPTIMIZER_VERBOSE 
#endif

#include <CGAL/Timer.h>
#include <CGAL/Origin.h>

#include <vector>
#include <list>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/format.hpp>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT,
          typename MoveFunction>
class Mesh_global_optimizer_2
{  
  // Types
  typedef CDT  Tr;
  typedef MoveFunction Mf;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Tr::Point            Point_2;
  typedef typename Tr::Face_handle      Face_handle;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::Vertex           Vertex;
  typedef typename Tr::Face_circulator  Face_circulator;

  typedef typename Gt::FT               FT;
  typedef typename Gt::Vector_2         Vector_2;
  
  typedef typename std::vector<Face_handle>                 Face_vector;
  typedef typename std::set<Vertex_handle>                  Vertex_set;
  typedef typename std::list<FT>                            FT_list;
  typedef typename std::pair<Vertex_handle,Point_2>         Move;

  typedef std::vector<Move>   Moves_vector;

  typedef typename MoveFunction::Sizing_field Sizing_field;

public:
  /**
   * Constructor
   */
  Mesh_global_optimizer_2(CDT& cdt,
                        const FT& freeze_ratio = 0.01,
                        const FT& convergence_ratio = 0.01,
                        const MoveFunction move_function = MoveFunction())
    : cdt_(cdt)
    , sq_freeze_ratio_(freeze_ratio * freeze_ratio)
    , convergence_ratio_(convergence_ratio)
    , move_function_(move_function)
    , sizing_field_(cdt)
  {}
  
  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }

  void operator()(const int nb_iterations)
  {
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
    CGAL::Timer timer;
    timer.start();
#endif
    // Fill set containing moving vertices
    Vertex_set moving_vertices;
    for(typename Tr::Finite_vertices_iterator
      vit = cdt_.finite_vertices_begin();
      vit != cdt_.finite_vertices_end();
      ++vit )
    {
      if(!cdt_.are_there_incident_constraints(vit))
        moving_vertices.insert(vit);
    }

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  double initial_vertices_nb = static_cast<double>(moving_vertices.size());
  double step_begin = timer.time();
  std::cerr << "Running " << Mf::name() << "-smoothing..." << std::endl;
  std::cerr << "(" << initial_vertices_nb << " vertices moving)" << std::endl;
#endif

    // Initialize big moves (stores the largest moves)
    big_moves_.clear();
    std::size_t big_moves_size = (std::max)(std::size_t(1),
                                            moving_vertices.size()/200);
    big_moves_.resize(big_moves_size, FT(0));

    // Iterate
    int i = -1;
    while ( ++i < nb_iterations && ! is_time_limit_reached() )
    {
      // Compute move for each vertex
      Moves_vector moves = compute_moves(moving_vertices);

      // Stop if convergence or time_limit is reached
      if ( check_convergence() || is_time_limit_reached() )
        break;

      // Update mesh with those moves
      update_mesh(moves, moving_vertices);

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
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
    
      step_begin = timer.time();
#endif
    }

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
    timer.stop();
    std::cerr << std::endl;
    if ( check_convergence() )
      std::cerr << "Convergence reached" << std::endl;

    std::cerr << "Total optimization time: " << timer.time()
              << "s" << std::endl << std::endl;
#endif
  }

private:
  /**
   * Returns moves for vertices of set \c moving_vertices
   */
  Moves_vector compute_moves(const Vertex_set& moving_vertices)
  {
    typename Gt::Construct_translated_point_2 translate =
      Gt().construct_translated_point_2_object();

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
      Vector_2 move = compute_move(*vit);
      if ( CGAL::NULL_VECTOR != move )
      {
        Point_2 new_position = translate((*vit)->point(), move);
        moves.push_back(std::make_pair(*vit,new_position));
      }

      // Stop if time_limit_ is reached
      if ( is_time_limit_reached() )
        break;
    }
    return moves;
  }

  /**
   * Returns the move for vertex \c v
   */
  Vector_2 compute_move(const Vertex_handle& v)
  {
    // Get move from move function
    Vector_2 move = move_function_(v, cdt_, sizing_field_);

    FT local_sq_size = min_sq_circumradius(v);
    if ( FT(0) == local_sq_size )
      return CGAL::NULL_VECTOR;

    FT local_move_sq_length = (move * move);

    // Move point only if displacement is big enough w.r.t local size
    if ( local_move_sq_length/local_sq_size < sq_freeze_ratio_ )
      return CGAL::NULL_VECTOR;

    // Update big moves
    update_big_moves(local_move_sq_length);

    return move;
  }

  /**
   * Returns the minimum cicumradius length of faces incident to \c v
   */
  FT min_sq_circumradius(const Vertex_handle& v) const
  {
    CGAL_assertion(!cdt_.is_infinite(v));

    Face_circulator face = cdt_.incident_faces(v);
    Face_circulator end = face;

    // Get first face sq_circumradius_length 
    // Initialize min
    FT min_sqr = std::numeric_limits<double>::max();
    // Find the minimum value
    do
    {
      if(face->is_in_domain())
        min_sqr = (std::min)(min_sqr, sq_circumradius(face));
      face++;
    }
    while(face != end);

    return min_sqr;
  }

  FT sq_circumradius(const Face_handle& f) const
  {
    Point_2 cc = cdt_.circumcenter(f);
    Point_2 p0 = f->vertex(0)->point();
    return (cc - p0) * (cc - p0);
  }

  /**
   * update big_moves_ vector with new_sq_move value
   */
  void update_big_moves(const FT& new_sq_move)
  {
    if ( new_sq_move > big_moves_.back() )
    {
      // Remove last value
      big_moves_.pop_back();

      // Insert value at the right place
      typename FT_list::iterator pos = std::find_if(
        big_moves_.begin(),
        big_moves_.end(),
        boost::lambda::_1 < new_sq_move );

      big_moves_.insert(pos, new_sq_move);
    }
  }

  bool is_time_limit_reached() const
  {
    return ( (time_limit() > 0) && (running_time_.time() > time_limit()) );      
  }

  bool check_convergence() const
  {
    namespace bl = boost::lambda;

    FT sum(0);
    for(typename FT_list::const_iterator it = big_moves_.begin();
        it != big_moves_.end();
        ++it)
      sum += CGAL::sqrt(*it);

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
    sum_moves_ = sum/big_moves_.size();
#endif

    return ( sum/FT(big_moves_.size()) < convergence_ratio_ );
  }

  void update_mesh(const Moves_vector& moves,
                   Vertex_set& moving_vertices)
  {
    // Facets that have to be updated
    std::set<Face_handle> outdated_faces;

    // Apply moves in triangulation
    for(typename Moves_vector::const_iterator it = moves.begin() ;
        it != moves.end() ;
        ++it )
    {
      const Vertex_handle& v = it->first;
      const Point_2& new_position = it->second;

      //cdt_.move(v, new_position);
      //function not available, see Constrained_triangulation_2
      cdt_.remove(v);
      cdt_.insert(new_position);
      //todo : insert new faces to outdated_faces

      if( is_time_limit_reached() )
        break;
    }

    // Update cdt
    moving_vertices.clear();
    reset_faces(outdated_faces.begin(), outdated_faces.end());
  }

  template<typename FaceIterator>
  void reset_faces(FaceIterator begin,
                   FaceIterator end)
  {
    ;//ToDo
    // update "inside" and "blind"
  }

private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  CDT& cdt_;
  FT sq_freeze_ratio_;
  FT convergence_ratio_;
  MoveFunction move_function_;
  Sizing_field sizing_field_;

  double time_limit_;
  CGAL::Timer running_time_;
  
  std::list<FT> big_moves_;
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  mutable FT sum_moves_;
#endif

};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_2_H
