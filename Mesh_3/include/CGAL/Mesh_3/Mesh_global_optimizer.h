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

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Real_timer.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>
#include <CGAL/Origin.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/Null_global_optimizer_visitor.h>
#include <CGAL/iterator.h>
#include <CGAL/tuple.h>

#include <CGAL/Mesh_3/Concurrent_mesher_config.h>

#ifdef CGAL_MESH_3_PROFILING
  #include <CGAL/Mesh_3/Profiling_tools.h>
#endif

#include <vector>
#include <list>
#include <limits>

#include <boost/type_traits/is_convertible.hpp>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/atomic.h>
# include <tbb/parallel_do.h>
# include <tbb/concurrent_vector.h>
#endif

namespace CGAL {

namespace Mesh_3 {


/************************************************
// Class Mesh_global_optimizer_base
// Two versions: sequential / parallel
************************************************/

// Sequential
template <typename Tr, typename Concurrency_tag>
class Mesh_global_optimizer_base
{
protected:
  typedef typename Tr::Geom_traits                          Gt;
  typedef typename Gt::FT                                   FT;
  typedef typename Tr::Lock_data_structure                  Lock_data_structure;

  typedef std::vector<cpp11::tuple<
    typename Tr::Vertex_handle, typename Tr::Weighted_point, FT> > Moves_vector;
  typedef unsigned int                                             Nb_frozen_points_type ;

  Mesh_global_optimizer_base(const Bbox_3 &, int)
    : big_moves_size_(0) {}

  void update_big_moves(const FT& new_sq_move)
  {
    if (big_moves_.size() < big_moves_size_ )
      big_moves_.insert(new_sq_move);
    else
    {
      FT smallest = *(big_moves_.begin());
      if( new_sq_move > smallest )
      {
        big_moves_.erase(big_moves_.begin());
        big_moves_.insert(new_sq_move);
      }
    }
  }

  void clear_big_moves()
  {
    big_moves_.clear();
  }

  Lock_data_structure *get_lock_data_structure() { return 0; }
  void unlock_all_elements() {}

protected:
  std::size_t big_moves_size_;
  std::multiset<FT> big_moves_;
};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template <typename Tr>
class Mesh_global_optimizer_base<Tr, Parallel_tag>
{
protected:
  typedef typename Tr::Geom_traits                          Gt;
  typedef typename Gt::FT                                   FT;
  typedef typename Tr::Lock_data_structure                  Lock_data_structure;
  typedef tbb::concurrent_vector<cpp11::tuple<
    typename Tr::Vertex_handle, typename Tr::Weighted_point, FT> >   Moves_vector;
  typedef tbb::atomic<unsigned int>                         Nb_frozen_points_type ;

  Mesh_global_optimizer_base(const Bbox_3 &bbox, int num_grid_cells_per_axis)
    : big_moves_size_(0)
    , m_lock_ds(bbox, num_grid_cells_per_axis)
  {
    big_moves_current_size_ = 0;
    big_moves_smallest_ = std::numeric_limits<FT>::max();
  }

  void update_big_moves(const FT& new_sq_move)
  {
    if (++big_moves_current_size_ <= big_moves_size_ )
    {
      tbb::mutex::scoped_lock lock(m_big_moves_mutex);
      typename std::multiset<FT>::const_iterator it = big_moves_.insert(new_sq_move);

      // New smallest move of all big moves?
      if (it == big_moves_.begin())
        big_moves_smallest_ = new_sq_move;
    }
    else
    {
      --big_moves_current_size_;

      if( new_sq_move > big_moves_smallest_ )
      {
        tbb::mutex::scoped_lock lock(m_big_moves_mutex);
        // Test it again since it may have been modified by another
        // thread in the meantime
        if( new_sq_move > big_moves_smallest_ )
        {
          big_moves_.erase(big_moves_.begin());
          typename std::multiset<FT>::const_iterator it = big_moves_.insert(new_sq_move);

          // New smallest move of all big moves?
          if (it == big_moves_.begin())
            big_moves_smallest_ = new_sq_move;
        }
      }
    }
  }

  void clear_big_moves()
  {
    big_moves_current_size_ = 0;
    big_moves_smallest_ = std::numeric_limits<FT>::max();
    big_moves_.clear();
  }

  Lock_data_structure *get_lock_data_structure()
  {
    return &m_lock_ds;
  }

  void unlock_all_elements()
  {
    m_lock_ds.unlock_all_points_locked_by_this_thread();
  }

public:

protected:
  tbb::atomic<std::size_t>  big_moves_current_size_;
  tbb::atomic<FT>           big_moves_smallest_;
  std::size_t               big_moves_size_;
  std::multiset<FT>         big_moves_;
  tbb::mutex                m_big_moves_mutex;

  /// Lock data structure
  Lock_data_structure m_lock_ds;
};
#endif // CGAL_LINKED_WITH_TBB




/************************************************
// Class Mesh_global_optimizer
************************************************/

template <typename C3T3,
          typename MeshDomain,
          typename MoveFunction,
          typename Visitor_ = Null_global_optimizer_visitor<C3T3> >
class Mesh_global_optimizer
: public Mesh_global_optimizer_base<typename C3T3::Triangulation, typename C3T3::Concurrency_tag>
{
  // Types
  typedef typename C3T3::Concurrency_tag Concurrency_tag;

  typedef Mesh_global_optimizer<C3T3, MeshDomain, MoveFunction, Visitor_> Self;
  typedef Mesh_global_optimizer_base<
    typename C3T3::Triangulation, typename C3T3::Concurrency_tag>         Base;

  using Base::get_lock_data_structure;
  using Base::big_moves_;
  using Base::big_moves_size_;

  typedef typename C3T3::Triangulation  Tr;
  typedef typename Tr::Geom_traits      Gt;

  typedef typename Tr::Bare_point       Bare_point;
  typedef typename Tr::Weighted_point   Weighted_point;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::Vertex           Vertex;

  typedef typename Gt::FT               FT;
  typedef typename Gt::Vector_3         Vector_3;

  typedef typename std::vector<Cell_handle>                 Cell_vector;
  typedef typename std::vector<Vertex_handle>               Vertex_vector;
  typedef typename std::set<Vertex_handle>                  Vertex_set;
  typedef typename Base::Moves_vector                       Moves_vector;
  typedef typename Base::Nb_frozen_points_type              Nb_frozen_points_type;

#ifdef CGAL_INTRUSIVE_LIST
  typedef Intrusive_list<Cell_handle>   Outdated_cell_set;
#else
  typedef std::set<Cell_handle>         Outdated_cell_set;
#endif //CGAL_INTRUSIVE_LIST

#ifdef CGAL_INTRUSIVE_LIST
  typedef Intrusive_list<Vertex_handle>  Moving_vertices_set;
#else
  typedef Vertex_set Moving_vertices_set;
#endif

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
                        const bool do_freeze,
                        const FT& convergence_ratio,
                        const MoveFunction move_function = MoveFunction());

  /**
   * Launch optimization process
   *
   * @param nb_interations maximum number of iterations
   */
  Mesh_optimization_return_code operator()(int nb_iterations,
                                           Visitor v = Visitor());

  /**
  * collects all vertices of the triangulation in moving_vertices
  * (even the frozen ones)
  */
  void collect_all_vertices(Moving_vertices_set& moving_vertices);

  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }

private:
  /**
   * Returns moves for vertices of set \c moving_vertices
   */
  Moves_vector compute_moves(Moving_vertices_set& moving_vertices)
{
  // Store new position of points which have to move
  Moves_vector moves;

  moves.reserve(moving_vertices.size());

  // reset worst_move list
  this->clear_big_moves();

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "Computing moves...";
  WallClockTimer t;
#endif


#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    tbb::concurrent_vector<Vertex_handle> vertices_not_moving_any_more;

    // Get move for each moving vertex
    tbb::parallel_do(
      moving_vertices.begin(), moving_vertices.end(),
      Compute_move<Self, Sizing_field, Moves_vector>(
        *this, sizing_field_, moves, do_freeze_, vertices_not_moving_any_more,
        tr_.geom_traits())
    );

    typename tbb::concurrent_vector<Vertex_handle>::const_iterator it
      = vertices_not_moving_any_more.begin();
    typename tbb::concurrent_vector<Vertex_handle>::const_iterator it_end
      = vertices_not_moving_any_more.end();
    for ( ; it != it_end ; ++it)
    {
      moving_vertices.erase(*it);
    }
  }
  // Sequential
  else
#endif // CGAL_LINKED_WITH_TBB
  {
    typename Gt::Construct_point_3 wp2p =
        tr_.geom_traits().construct_point_3_object();
    typename Gt::Construct_weighted_point_3 p2wp =
        tr_.geom_traits().construct_weighted_point_3_object();
    typename Gt::Construct_translated_point_3 translate =
        tr_.geom_traits().construct_translated_point_3_object();

    // Get move for each moving vertex
    typename Moving_vertices_set::iterator vit = moving_vertices.begin();
    for ( ; vit != moving_vertices.end() ; )
    {
      Vertex_handle oldv = *vit;
      ++vit;
      Vector_3 move = compute_move(oldv);

      if ( CGAL::NULL_VECTOR != move )
      {
        Bare_point new_position = translate(wp2p(oldv->point()),move);
        FT size = (Sizing_field::is_vertex_update_needed ?
          sizing_field_(new_position, oldv) : 0);
        moves.push_back(cpp11::make_tuple(oldv,p2wp(new_position),size));
      }
      else // CGAL::NULL_VECTOR == move
      {
        if(do_freeze_)
          moving_vertices.erase(oldv); // TODO: if non-intrusive,
                                       // we can optimize since we have the iterator,
                                       // don't forget to do "vit = mv.erase(vit)" instead ++vit
      }

      // Stop if time_limit_ is reached
      if ( is_time_limit_reached() )
        break;
    }
  }

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "done in " << t.elapsed() << " seconds." << std::endl;
#endif

  return moves;
}


  /**
   * Returns the move for vertex \c v
   * warning : this function should be called only on moving vertices
   *           even for frozen vertices, it could return a non-zero vector
   */
  Vector_3 compute_move(const Vertex_handle& v);

  /**
   * Updates mesh using moves of \c moves vector. Updates moving_vertices with
   * the new set of moving vertices after the move.
   */
  void update_mesh(const Moves_vector& moves,
                   Moving_vertices_set& moving_vertices,
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
  FT min_circumradius_sq_length(const Vertex_handle& v, const Cell_vector& incident_cells) const;

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

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for compute_moves function
  template <typename MGO, typename Sizing_field_, typename Moves_vector_>
  class Compute_move
  {
    typedef tbb::concurrent_vector<Vertex_handle> Vertex_conc_vector;

    MGO                  & m_mgo;
    const Sizing_field_  & m_sizing_field;
    Moves_vector_        & m_moves;
    bool                   m_do_freeze;
    Vertex_conc_vector   & m_vertices_not_moving_any_more;
    const Gt             & m_gt;

  public:
    // Constructor
    Compute_move(MGO &mgo, 
                 const Sizing_field_ &sizing_field,
                 Moves_vector_ &moves,
                 bool do_freeze,
                 Vertex_conc_vector &vertices_not_moving_any_more,
                 const Gt &gt)
    : m_mgo(mgo),
      m_sizing_field(sizing_field),
      m_moves(moves),
      m_do_freeze(do_freeze),
      m_vertices_not_moving_any_more(vertices_not_moving_any_more),
      m_gt(gt)
    {}

    // Constructor
    Compute_move(const Compute_move &cm)
    : m_mgo(cm.m_mgo), 
      m_sizing_field(cm.m_sizing_field), 
      m_moves(cm.m_moves),
      m_do_freeze(cm.m_do_freeze),
      m_vertices_not_moving_any_more(cm.m_vertices_not_moving_any_more),
      m_gt(cm.m_gt)
    {}

    // operator()
    void operator()(const Vertex_handle& oldv) const
    {
      typename Gt::Construct_point_3 wp2p =
          m_gt.construct_point_3_object();
      typename Gt::Construct_weighted_point_3 p2wp =
          m_gt.construct_weighted_point_3_object();
      typename Gt::Construct_translated_point_3 translate =
          m_gt.construct_translated_point_3_object();

      Vector_3 move = m_mgo.compute_move(oldv);
      if ( CGAL::NULL_VECTOR != move )
      {
        Bare_point new_position = translate(wp2p(oldv->point()), move);
        FT size = (MGO::Sizing_field::is_vertex_update_needed ?
          m_sizing_field(new_position, oldv) : 0);
        // typedef Triangulation_helpers<typename C3T3::Triangulation> Th;
        //if( !Th().inside_protecting_balls(tr_, oldv, new_position))
        //note : this is not happening for Lloyd and ODT so it's commented
        //       maybe for a new global optimizer it should be de-commented
        m_moves.push_back(cpp11::make_tuple(oldv, p2wp(new_position), size));
      }
      else // CGAL::NULL_VECTOR == move
      {
        if(m_do_freeze)
        {
          m_vertices_not_moving_any_more.push_back(oldv);
        }
      }

      if ( m_mgo.is_time_limit_reached() )
        tbb::task::self().cancel_group_execution();
    }
  };

  // Functor for fill_sizing_field function
  template <typename MGO, typename Local_list_>
  class Compute_sizing_field_value
  {
    MGO                  & m_mgo;
    const Gt             & m_gt;
    Local_list_          & m_local_lists;

  public:
    // Constructor
    Compute_sizing_field_value(MGO &mgo, 
                               const Gt &gt,
                               Local_list_ &local_lists)
    : m_mgo(mgo),
      m_gt(gt),
      m_local_lists(local_lists)
    {}

    // Constructor
    Compute_sizing_field_value(const Compute_sizing_field_value &csfv)
    : m_mgo(csfv.m_mgo), 
      m_gt(csfv.m_gt),
      m_local_lists(csfv.m_local_lists)
    {}

    // operator()
    void operator()(Vertex& v) const
    {
      typename Gt::Construct_point_3 wp2p =
          m_gt.construct_point_3_object();

      Vertex_handle vh 
        = Tr::Triangulation_data_structure::Vertex_range::s_iterator_to(v);
      m_local_lists.local().push_back(
          std::make_pair(wp2p(v.point()), m_mgo.average_circumradius_length(vh)));
    }
  };

  // Functor for update_mesh function
  template <typename MGO, typename Helper, typename Tr_, typename Moves_vector_,
            typename Moving_vertices_set_, typename Outdated_cell_set_>
  class Move_vertex
  {
    MGO                                  & m_mgo;
    const Helper                         & m_helper;
    const Moves_vector_                  & m_moves;
    Moving_vertices_set_                 & m_moving_vertices;
    Outdated_cell_set_                   & m_outdated_cells;
  
    typedef typename Tr_::Bare_point    Bare_point;
    typedef typename Tr_::Vertex_handle Vertex_handle;

  public:
    // Constructor
    Move_vertex(MGO &mgo, const Helper &helper, const Moves_vector_ &moves,
                Moving_vertices_set_ &moving_vertices,
                Outdated_cell_set_ &outdated_cells)
    : m_mgo(mgo), m_helper(helper), m_moves(moves), 
      m_moving_vertices(moving_vertices), m_outdated_cells(outdated_cells)
    {}

    // Constructor
    Move_vertex(const Move_vertex &mv)
    : m_mgo(mv.m_mgo), m_helper(mv.m_helper), m_moves(mv.m_moves),
      m_moving_vertices(mv.m_moving_vertices),
      m_outdated_cells(mv.m_outdated_cells)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
      {
        const Vertex_handle& v = cpp11::get<0>(m_moves[i]);
        const Weighted_point& new_position = cpp11::get<1>(m_moves[i]);
        // Get size at new position
        if ( MGO::Sizing_field::is_vertex_update_needed )
        {
          //FT size = sizing_field_(new_position,v);
          FT size = cpp11::get<2>(m_moves[i]);

          // Move point
          bool could_lock_zone;
          Vertex_handle new_v = m_helper.move_point(
            v, new_position, m_outdated_cells, m_moving_vertices, &could_lock_zone);
          while (could_lock_zone == false)
          {
            new_v = m_helper.move_point(
              v, new_position, m_outdated_cells, m_moving_vertices, &could_lock_zone);
          }

          // Restore size in meshing_info data
          new_v->set_meshing_info(size);
        }
        else // Move point
        {
          bool could_lock_zone;
          do {
            m_helper.move_point(
              v, new_position, m_outdated_cells, m_moving_vertices, &could_lock_zone);
          } while (!could_lock_zone);
        }

        m_mgo.unlock_all_elements();

        // Stop if time_limit_ is reached, here we can't return without rebuilding
        // restricted delaunay
        if ( m_mgo.is_time_limit_reached() )
        {
          tbb::task::self().cancel_group_execution();
          break;
        }
      }
    }
  };
#endif // CGAL_LINKED_WITH_TBB

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
  CGAL::Real_timer running_time_;

  bool do_freeze_;
  mutable Nb_frozen_points_type nb_frozen_points_;

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  mutable FT sum_moves_;
#endif
};


template <typename C3T3, typename Md, typename Mf, typename V_>
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
Mesh_global_optimizer(C3T3& c3t3,
                      const Md& domain,
                      const FT& freeze_ratio,
                      const bool do_freeze,
                      const FT& convergence_ratio,
                      const Mf move_function)
: Base(c3t3.bbox(),
       Concurrent_mesher_config::get().locking_grid_num_cells_per_axis)
, c3t3_(c3t3)
, tr_(c3t3_.triangulation())
, domain_(domain)
, sq_freeze_ratio_(freeze_ratio*freeze_ratio)
, convergence_ratio_(convergence_ratio)
, helper_(c3t3_,domain_, get_lock_data_structure())
, move_function_(move_function)
, sizing_field_(c3t3.triangulation())
, time_limit_(-1)
, running_time_()

, do_freeze_(do_freeze)

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
, sum_moves_(0)
#endif // CGAL_MESH_3_OPTIMIZER_VERBOSE
{
  nb_frozen_points_ = 0; // We put it here in case it's an "atomic"

  // If we're multi-thread
  tr_.set_lock_data_structure(get_lock_data_structure());

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  std::cerr << "Fill sizing field...";
  CGAL::Real_timer timer;
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

#ifdef CGAL_MESH_3_PROFILING
  WallClockTimer t;
#endif

  // Fill set containing moving vertices
  // first, take them all
#ifdef  CGAL_CONSTRUCT_INTRUSIVE_LIST_RANGE_CONSTRUCTOR
  typedef CGAL::Prevent_deref<typename Tr::Finite_vertices_iterator> It;
  Moving_vertices_set moving_vertices(It(tr_.finite_vertices_begin()),
                                      It(tr_.finite_vertices_end()));
#else
  Moving_vertices_set moving_vertices;
  collect_all_vertices(moving_vertices);
#endif

  std::size_t initial_vertices_nb = moving_vertices.size();
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  double step_begin = running_time_.time();
  std::cerr << "Running " << Mf::name() << "-smoothing ("
    << initial_vertices_nb << " vertices)" << std::endl;
#endif //CGAL_MESH_3_OPTIMIZER_VERBOSE

  // Initialize big moves (stores the largest moves)
  this->clear_big_moves();
  big_moves_size_ =
    (std::max)(std::size_t(1), std::size_t(moving_vertices.size()/500));

  std::size_t nb_vertices_moved = -1;
  bool convergence_stop = false;

  // Iterate
  int i = -1;
  while ( ++i < nb_iterations && ! is_time_limit_reached() )
  {
    if(!do_freeze_)
      nb_frozen_points_ = 0;
    else
      nb_vertices_moved = moving_vertices.size();

    // Compute move for each vertex
    Moves_vector moves = compute_moves(moving_vertices);
    visitor.after_compute_moves();

    //Pb with Freeze : sometimes a few vertices continue moving indefinitely
    //if the nb of moving vertices is < 1% of total nb AND does not decrease
    if(do_freeze_
      && double(nb_vertices_moved) < 0.005 * double(initial_vertices_nb)
      && nb_vertices_moved == moving_vertices.size())
    {
      // we should stop because we are
      // probably entering an infinite instable loop
      convergence_stop = true;
      break;
    }

    // Stop if time_limit is reached
    if ( is_time_limit_reached() )
      break;

    // Update mesh with those moves
    update_mesh(moves, moving_vertices, visitor);
    visitor.end_of_iteration(i);

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
    std::size_t moving_vertices_size = moving_vertices.size();
    std::cerr << boost::format("\r             \r"
      "end iter.%1%, %2% / %3%, last step:%4$.2fs, step avg:%5$.2fs, avg big moves:%6$.3f ")
    % (i+1)
    % moving_vertices_size
    % initial_vertices_nb
    % (running_time_.time() - step_begin)
    % (running_time_.time() / (i+1))
    % sum_moves_;
    step_begin = running_time_.time();
#endif

    if (do_freeze_ && nb_frozen_points_ == initial_vertices_nb )
      break;

    if(check_convergence())
      break;
  }
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  std::cerr << std::endl;
#endif

  running_time_.stop();

#ifdef CGAL_MESH_3_PROFILING
  double optim_time = t.elapsed();
# ifdef CGAL_MESH_3_EXPORT_PERFORMANCE_DATA
  CGAL_MESH_3_SET_PERFORMANCE_DATA(std::string(Mf::name()) + "_optim_time", optim_time);
# endif
#endif

#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  if ( do_freeze_ && nb_frozen_points_ == initial_vertices_nb )
    std::cerr << "All vertices frozen" << std::endl;
  else if ( do_freeze_ && convergence_stop )
    std::cerr << "Can't improve anymore" << std::endl;
  else if ( is_time_limit_reached() )
    std::cerr << "Time limit reached" << std::endl;
  else if ( check_convergence() )
    std::cerr << "Convergence reached" << std::endl;
  else if( i >= nb_iterations )
    std::cerr << "Max iteration number reached" << std::endl;

  std::cerr << "Total optimization time: " << running_time_.time()
            << "s" << std::endl << std::endl;
#endif

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << std::endl << "Total optimization 'wall-clock' time: "
            << optim_time << "s" << std::endl;
#endif

  if ( do_freeze_ && nb_frozen_points_ == initial_vertices_nb )
    return ALL_VERTICES_FROZEN;

  else if ( do_freeze_ && convergence_stop )
    return CANT_IMPROVE_ANYMORE;

  else if ( is_time_limit_reached() )
    return TIME_LIMIT_REACHED;

  else if ( check_convergence() )
    return CONVERGENCE_REACHED;

  return MAX_ITERATION_NUMBER_REACHED;
}

template <typename C3T3, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
collect_all_vertices(Moving_vertices_set& moving_vertices)
{
  typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
  for(; vit != tr_.finite_vertices_end(); ++vit )
    moving_vertices.insert(vit);
}

template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::Vector_3
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
compute_move(const Vertex_handle& v)
{
  typename Gt::Compute_squared_length_3 sq_length =
      tr_.geom_traits().compute_squared_length_3_object();
  typename Gt::Construct_vector_3 vector =
      tr_.geom_traits().construct_vector_3_object();
  typename Gt::Construct_translated_point_3 translate =
      tr_.geom_traits().construct_translated_point_3_object();
  typename Gt::Construct_point_3 wp2p =
      tr_.geom_traits().construct_point_3_object();

  Cell_vector incident_cells;
  incident_cells.reserve(64);
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    tr_.incident_cells_threadsafe(v, std::back_inserter(incident_cells));
  }
  else
#endif //CGAL_LINKED_WITH_TBB
  {
    tr_.incident_cells(v, std::back_inserter(incident_cells));
  }

  // Get move from move function
  Vector_3 move = move_function_(v, incident_cells, c3t3_, sizing_field_);

  // Project surface vertex
  if ( c3t3_.in_dimension(v) == 2 )
  {
    Bare_point new_position = translate(wp2p(v->point()),move);
    move = vector(wp2p(v->point()), helper_.project_on_surface(new_position, v));
  }

  FT local_sq_size = min_circumradius_sq_length(v, incident_cells);
  if ( FT(0) == local_sq_size )
    return CGAL::NULL_VECTOR;

  FT local_move_sq_ratio = sq_length(move) / local_sq_size;

  // Move point only if displacement is big enough w.r.t local size
  if ( local_move_sq_ratio < sq_freeze_ratio_ )
  {
    nb_frozen_points_++;
    return CGAL::NULL_VECTOR;
  }

  // Update big moves
  this->update_big_moves(local_move_sq_ratio);

  return move;
}


template <typename C3T3, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
update_mesh(const Moves_vector& moves,
            Moving_vertices_set& moving_vertices,
            Visitor& visitor)
{
  // Cells which have to be updated
  Outdated_cell_set outdated_cells;

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "Moving vertices...";
  WallClockTimer t;
#endif

#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    // Apply moves in triangulation
    tbb::parallel_for(tbb::blocked_range<size_t>(0, moves.size()),
      Move_vertex<
        Self, C3T3_helpers, Tr, Moves_vector,
        Moving_vertices_set, Outdated_cell_set>(
          *this, helper_, moves, moving_vertices, outdated_cells)
    );
  }
  // Sequential
  else
#endif // CGAL_LINKED_WITH_TBB
  {
    // Apply moves in triangulation
    for ( typename Moves_vector::const_iterator it = moves.begin() ;
         it != moves.end() ;
         ++it )
    {
      const Vertex_handle& v = cpp11::get<0>(*it);
      const Weighted_point& new_position = cpp11::get<1>(*it);
      // Get size at new position
      if ( Sizing_field::is_vertex_update_needed )
      {
        //FT size = sizing_field_(new_position,v);
        FT size = cpp11::get<2>(*it);

        // Move point
        Vertex_handle new_v = helper_.move_point(v, new_position, outdated_cells, moving_vertices);

        // Restore size in meshing_info data
        new_v->set_meshing_info(size);
      }
      else // Move point
      {
        helper_.move_point(v, new_position, outdated_cells, moving_vertices);
      }

      // Stop if time_limit_ is reached, here we can't return without rebuilding
      // restricted delaunay
      if ( is_time_limit_reached() )
        break;
    }
  }


  visitor.after_move_points();

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "done in " << t.elapsed() << " seconds." << std::endl;
#endif


#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "Updating C3T3 (rebuilding restricted Delaunay)...";
  t.reset();
#endif

  // Update c3t3
#ifdef CGAL_INTRUSIVE_LIST
  // AF:  rebuild_restricted_delaunay does more cell insertion/removal
  //      which clashes with our inplace list
  //      That's why we hand it in, and call clear() when we no longer need it
  helper_.rebuild_restricted_delaunay(outdated_cells,
                                      moving_vertices);
#else
  helper_.rebuild_restricted_delaunay(outdated_cells.begin(),
                                      outdated_cells.end(),
                                      moving_vertices);
#endif

  visitor.after_rebuild_restricted_delaunay();

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "Updating C3T3 done in " << t.elapsed() << " seconds." << std::endl;
#endif

}


template <typename C3T3, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
fill_sizing_field()
{
  std::map<Bare_point,FT> value_map;

#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    typedef tbb::enumerable_thread_specific<
      std::vector< std::pair<Bare_point, FT> > > Local_list;
    Local_list local_lists;

    tbb::parallel_do(
      tr_.finite_vertices_begin(), tr_.finite_vertices_end(),
      Compute_sizing_field_value<Self, Local_list>(*this, tr_.geom_traits(), local_lists)
    );

    for(typename Local_list::iterator it_list = local_lists.begin() ;
          it_list != local_lists.end() ;
          ++it_list )
    {
      value_map.insert(it_list->begin(), it_list->end());
    }
  }
  else
#endif //CGAL_LINKED_WITH_TBB
  {
    // Fill map with local size
    for(typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
        vit != tr_.finite_vertices_end();
        ++vit)
    {
      typename Gt::Construct_point_3 wp2p =
          tr_.geom_traits().construct_point_3_object();

      value_map.insert(std::make_pair(wp2p(vit->point()),
                                      average_circumradius_length(vit)));
    }
  }

  // fill sizing field
  sizing_field_.fill(value_map);
}


template <typename C3T3, typename Md, typename Mf, typename V_>
bool
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
check_convergence() const
{
  FT sum(0);
  for( typename std::multiset<FT>::const_iterator
       it = big_moves_.begin(), end = big_moves_.end() ; it != end ; ++it )
  {
    sum += CGAL::sqrt(*it);
  }

  FT average_move = sum/FT(big_moves_size_);/*even if set is not full, divide*/
       /*by max size so that if only 1 point moves, it goes to 0*/
#ifdef CGAL_MESH_3_OPTIMIZER_VERBOSE
  sum_moves_ = average_move;
#endif

  return ( average_move < convergence_ratio_ );
}


template <typename C3T3, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C3T3,Md,Mf,V_>::FT
Mesh_global_optimizer<C3T3,Md,Mf,V_>::
average_circumradius_length(const Vertex_handle& v) const
{
  Cell_vector incident_cells;
  incident_cells.reserve(64);
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    tr_.incident_cells_threadsafe(v, std::back_inserter(incident_cells));
  }
  else
#endif //CGAL_LINKED_WITH_TBB
  {
    tr_.incident_cells(v, std::back_inserter(incident_cells));
  }

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
min_circumradius_sq_length(const Vertex_handle& v, const Cell_vector& incident_cells) const
{
  // Get first cell sq_circumradius_length
  typename Cell_vector::const_iterator cit = incident_cells.begin();
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
    tr_.geom_traits().compute_squared_distance_3_object();
  typename Gt::Construct_point_3 wp2p =
    tr_.geom_traits().construct_point_3_object();

  const Bare_point circumcenter = tr_.dual(cell);
  return ( sq_distance(wp2p(v->point()), circumcenter) );
}

} // end namespace Mesh_3

} //namespace CGAL

#endif // CGAL_MESH_3_MESH_GLOBAL_OPTIMIZER_H
