// Copyright (c) 2009-2013 INRIA Sophia-Antipolis (France),
//               2014-2015 GeometryFactory (France).
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
// Author(s) : Jane Tournois, Raul Gallegos, Pierre Alliez, Stéphane Tayeb
//

#ifndef CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_2_H
#define CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_2_H

#ifdef CGAL_MESH_2_VERBOSE
  #define CGAL_MESH_2_OPTIMIZER_VERBOSE 
#endif

#include <CGAL/Timer.h>
#include <CGAL/Origin.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_voronoi_diagram_2.h>

#include <vector>
#include <list>
#include <algorithm>
#include <iterator>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

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
                        const FT& convergence_ratio = 0., //no criterion
                        const FT& freeze_ratio = 0., //no criterion
                        const MoveFunction move_function = MoveFunction())
    : cdt_(cdt)
    , sq_freeze_ratio_(freeze_ratio * freeze_ratio)
    , convergence_ratio_(convergence_ratio)
    , move_function_(move_function)
    , sizing_field_(cdt)
    , seeds_()
    , seeds_mark_(false)
  {
  }
  
  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }

  /** The value type of \a InputIterator should be \c Point, and represents
      seeds.
  */
  template<typename InputIterator>
  void set_seeds(InputIterator b,
                 InputIterator e,
                 const bool mark = false)
  {
    if(b != e)
    {
      seeds_.clear();
      std::copy(b, e, std::back_inserter(seeds_));
      seeds_mark_ = mark;
    }
  }

  Mesh_optimization_return_code operator()(const int nb_iterations)
  {
    running_time_.reset();
    running_time_.start();

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

  double initial_vertices_nb = static_cast<double>(moving_vertices.size());
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  double step_begin = running_time_.time();
  std::cerr << "Running " << Mf::name() << "-smoothing..." << std::endl;
  std::cerr << "(" << initial_vertices_nb << " vertices moving)" << std::endl;
#endif

    // Initialize big moves (stores the largest moves)
    big_moves_.clear();
    std::size_t big_moves_size = (std::max)(std::size_t(1),
                                            moving_vertices.size()/100);
    big_moves_.resize(big_moves_size, FT(0));

    std::size_t nb_vertices_moved = -1;
    bool convergence_stop = false;

    // Iterate
    int i = -1;
    while ( ++i < nb_iterations && ! is_time_limit_reached() )
    {
      this->before_move();

      // Compute move for each vertex
      Moves_vector moves = compute_moves(moving_vertices);

      //Pb with Freeze : sometimes a few vertices continue moving indefinitely
      //if the nb of moving vertices is < 2% of total nb AND does not decrease
      if(sq_freeze_ratio_ > 0.
        && nb_vertices_moved < 0.01 * initial_vertices_nb
        && nb_vertices_moved == moving_vertices.size())
      { 
        // we should stop because we are 
        // probably entering an infinite instable loop
        convergence_stop = true;
        break;
      }

      // Stop if convergence or time_limit is reached
      if ( check_convergence() || is_time_limit_reached() )
        break;

      // Update mesh with those moves
      update_mesh(moves);
      nb_vertices_moved = moving_vertices.size();

      this->after_move();
 
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
      double time = running_time_.time();
      double moving_vertices_size = static_cast<double>(moving_vertices.size());
      std::cerr << boost::format("\r             \r"
        "end iteration %1% (%2%%% frozen), %3% / %4%, last step:%5$.2fs, step avg:%6$.2fs, avg large move:%7$.3f          ")
      % (i+1)
      % ((1. - moving_vertices_size/initial_vertices_nb)*100.)
      % moving_vertices_size
      % initial_vertices_nb
      % (time - step_begin)
      % (time / (i+1))
      % sum_moves_;
      step_begin = running_time_.time();
#endif
    }

    this->after_all_moves();
    running_time_.stop();

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
    if(sq_freeze_ratio_ > 0. && moving_vertices.empty())
      std::cerr << "All vertices frozen" << std::endl;
    else if(sq_freeze_ratio_ > 0. && convergence_stop)
      std::cerr << "Can't improve anymore" << std::endl;
    else if ( is_time_limit_reached() )
      std::cerr << "Time limit reached" << std::endl;
    else if ( check_convergence() )
      std::cerr << "Convergence reached" << std::endl;
    else if ( i >= nb_iterations )
      std::cerr << "Max iteration number reached" << std::endl;

    std::cerr << "Total optimization time: " << running_time_.time()
              << "s" << std::endl << std::endl;
#endif

    if( sq_freeze_ratio_ > 0. && moving_vertices.empty() )
      return ALL_VERTICES_FROZEN;
    else if( sq_freeze_ratio_ > 0. && convergence_stop )
      return CANT_IMPROVE_ANYMORE;
    else if( is_time_limit_reached() )
      return TIME_LIMIT_REACHED;
    else if( check_convergence() )
      return CONVERGENCE_REACHED;
    else if( i >= nb_iterations )
      return MAX_ITERATION_NUMBER_REACHED;
    else
      return MESH_OPTIMIZATION_UNKNOWN_ERROR;
  }

private:
  /**
   * Returns moves for vertices of set \c moving_vertices
   */
  Moves_vector compute_moves(Vertex_set& moving_vertices)
  {
    typename Gt::Construct_translated_point_2 translate =
      Gt().construct_translated_point_2_object();

    // Store new location of points which have to move
    Moves_vector moves;
    moves.reserve(moving_vertices.size());

    // reset worst_move list
    std::fill(big_moves_.begin(), big_moves_.end(), FT(0));

    // Get move for each moving vertex
    for ( typename Vertex_set::const_iterator vit = moving_vertices.begin() ;
      vit != moving_vertices.end() ; )
    {
      Vertex_handle oldv = *vit;
      Vector_2 move = compute_move(oldv);
      ++vit;

      if ( CGAL::NULL_VECTOR != move )
      {
        Point_2 new_position = translate(oldv->point(), move);
        moves.push_back(std::make_pair(oldv, new_position));
      }
      else if(sq_freeze_ratio_ > 0.) //freezing ON
        moving_vertices.erase(oldv);

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

    FT local_move_sq_ratio = (move * move) / local_sq_size;

    // Move point only if displacement is big enough w.r.t local size
    if ( local_move_sq_ratio < sq_freeze_ratio_ )
      return CGAL::NULL_VECTOR;

    // Update big moves
    update_big_moves(local_move_sq_ratio);

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
    FT min_sqr = (std::numeric_limits<double>::max)();
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
    return (time_limit() > 0)
        && (running_time_.time() > time_limit());
  }

  bool check_convergence() const
  {
    namespace bl = boost::lambda;

    FT sum(0);
    for(typename FT_list::const_iterator it = big_moves_.begin();
        it != big_moves_.end();
        ++it)
      sum += CGAL::sqrt(*it);

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
    sum_moves_ = sum/big_moves_.size();
#endif

    return ( sum/FT(big_moves_.size()) < convergence_ratio_ );
  }

  void update_mesh(const Moves_vector& moves)
  {
    // Apply moves in triangulation
    for(typename Moves_vector::const_iterator it = moves.begin() ;
        it != moves.end() ;
        ++it )
    {
      const Vertex_handle& v = it->first;
      const Point_2& new_position = it->second;
      const FT size = v->sizing_info();

      //cdt_.move(v, new_position);
      //function not available, see Constrained_triangulation_2
      cdt_.remove(v);
      Vertex_handle new_v = cdt_.insert(new_position);

      new_v->set_sizing_info(size);

      if( is_time_limit_reached() )
        break;
    }
  }

  void before_move()
  {
    update_blind_faces();
  }

  void after_move()
  {
    //update inside/outside tags
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
    CGAL::Delaunay_mesher_2<CDT, Criteria>::mark_facets(cdt_,
      seeds_.begin(),
      seeds_.end(),
      seeds_mark_/*faces that are not in domain are tagged false*/);
    //Connected components of seeds are marked with the value of
    //  \a mark. Other components are marked with \c !mark. The connected
    //  component of infinite faces is always marked with \c false.
  }

  void after_all_moves()
  {
    update_blind_faces();
  }

  void update_blind_faces()
  {
    //update blindness
    CGAL::Constrained_voronoi_diagram_2<CDT> cvd(cdt_);
    cvd.tag_faces_blind();
  }

public:

  void output_angles_histogram(std::ostream& os)
  {
    double min_angle = 180.;
    double max_angle = 0.;

    //fill histogram
    std::vector<int> histo(180,0);
    for(typename CDT::Finite_faces_iterator fit = cdt_.finite_faces_begin();
        fit != cdt_.finite_faces_end();
        ++fit)
    {
      if(!fit->is_in_domain())
        continue;
      typename CDT::Triangle tr = cdt_.triangle(fit);
      for(int i = 0; i < 3; ++i)
      {
        Vector_2 v1(tr[i], tr[(i+1)%3]);
        Vector_2 v2(tr[i], tr[(i+2)%3]);
        v1 = v1 / CGAL::sqrt(v1.squared_length());
        v2 = v2 / CGAL::sqrt(v2.squared_length());

        using namespace boost::math::constants;
        double angle = std::acos(v1 * v2) * 180. / pi<double>();
        histo[ (std::min)(179, int(angle)) ]++;

        min_angle = (std::min)(min_angle, angle);
        max_angle = (std::max)(max_angle, angle);
      }
    }

    //output histogram to os
    os << "# Min angle = " << min_angle << std::endl;
    os << "# Max angle = " << max_angle << std::endl;
    os << std::endl << "## Histogram ##" << std::endl;
    for(unsigned int i = 0; i < 180; ++i)
      os << i << "\t" << histo[i] << std::endl;
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
  std::list<Point_2> seeds_;
  bool seeds_mark_;

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
