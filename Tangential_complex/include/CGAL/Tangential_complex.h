// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin


#ifndef TANGENTIAL_COMPLEX_H
#define TANGENTIAL_COMPLEX_H

#include <CGAL/Tangential_complex/config.h>
#include <CGAL/Tangential_complex/Simplicial_complex.h>
#include <CGAL/Tangential_complex/utilities.h>
#include <CGAL/Tangential_complex/Point_cloud.h>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Dimension.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/point_generators_d.h>

# include <CGAL/Mesh_3/Profiling_tools.h>

#include <CGAL/IO/Triangulation_off_ostream.h> // CJTODO TEMP

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <boost/optional.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <set>
#include <utility>
#include <sstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <functional>
#include <iterator>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/parallel_for.h>
# include <tbb/combinable.h>
# include <tbb/mutex.h>
#endif

//#define CGAL_TC_EXPORT_NORMALS // Only for 3D surfaces (k=2, d=3)
//#define CGAL_ALPHA_TC
const double ALPHA = 0.3;
#ifdef CGAL_LINKED_WITH_TBB
tbb::atomic<unsigned int> ttt_star; // CJTODO TEMP
tbb::atomic<unsigned int> ttt_intersect;
#endif

//CJTODO: debug
//#define CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_SPHERE_3
//#define CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_TORUS_D

namespace CGAL {

using namespace Tangential_complex_;

enum Fix_inconsistencies_status { 
  TC_FIXED = 0, TIME_LIMIT_REACHED, FIX_NOT_PERFORMED };

class Vertex_data
{
public:
  Vertex_data(std::size_t data = std::numeric_limits<std::size_t>::max())
    : m_data(data)
  {}
  operator std::size_t() { return m_data; }
  operator std::size_t() const { return m_data; }

private:
  std::size_t m_data;
};

/// The class Tangential_complex represents a tangential complex
template <
  typename Kernel, // ambiant dimension
  typename DimensionTag, // intrinsic dimension
  typename Concurrency_tag = CGAL::Parallel_tag,
  typename Tr = Regular_triangulation
  <
    Regular_triangulation_euclidean_traits<
      Epick_d<DimensionTag> >,

    Triangulation_data_structure
    <
      typename Regular_triangulation_euclidean_traits<
        Epick_d<DimensionTag> >::Dimension,
      Triangulation_vertex<Regular_triangulation_euclidean_traits<
        Epick_d<DimensionTag> >, Vertex_data >,
      Triangulation_full_cell<Regular_triangulation_euclidean_traits<
        Epick_d<DimensionTag> > >
    >
  >
>
class Tangential_complex
{
  typedef typename Kernel::FT                         FT;
  typedef typename Kernel::Point_d                    Point;
  typedef typename Kernel::Weighted_point_d           Weighted_point;
  typedef typename Kernel::Vector_d                   Vector;

  typedef Tr                                          Triangulation;
  typedef typename Triangulation::Geom_traits         Tr_traits;
  typedef typename Triangulation::Weighted_point      Tr_point;
  typedef typename Triangulation::Bare_point          Tr_bare_point;
  typedef typename Triangulation::Vertex_handle       Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle    Tr_full_cell_handle;
  typedef typename Tr_traits::Vector_d                Tr_vector;

  typedef Basis<Kernel>                               Tangent_space_basis;
  typedef Basis<Kernel>                               Orthogonal_space_basis;

  typedef std::vector<Point>                          Points;
#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
  typedef tbb::mutex                                  Mutex_for_perturb;
  typedef Vector                                      Translation_for_perturb;
  typedef std::vector<Atomic_wrapper<FT> >            Weights;
 #ifdef CGAL_TC_PERTURB_WEIGHT
  typedef std::vector<Atomic_wrapper<FT> >            Weights_memory;
 #endif
#else
  typedef Vector                                      Translation_for_perturb;
  typedef std::vector<FT>                             Weights;
 #ifdef CGAL_TC_PERTURB_WEIGHT
  typedef std::vector<FT>                             Weights_memory;
 #endif
#endif
  typedef std::vector<Translation_for_perturb>        Translations_for_perturb;

  typedef Point_cloud_data_structure<Kernel, Points>  Points_ds;
  typedef typename Points_ds::KNS_range               KNS_range;
  typedef typename Points_ds::KNS_iterator            KNS_iterator;
  typedef typename Points_ds::INS_range               INS_range;
  typedef typename Points_ds::INS_iterator            INS_iterator;

  // Store a local triangulation and a handle to its center vertex
  struct Tr_and_VH
  {
  public:
    Tr_and_VH()
      : m_tr(NULL) {}
    Tr_and_VH(int dim)
      : m_tr(new Triangulation(dim)) {}

    ~Tr_and_VH() { destroy_triangulation(); }

    Triangulation & construct_triangulation(int dim)
    {
      delete m_tr;
      m_tr = new Triangulation(dim);
      return tr();
    }

    void destroy_triangulation()
    {
      delete m_tr;
      m_tr = NULL;
    }

    Triangulation &      tr()       { return *m_tr; }
    Triangulation const& tr() const { return *m_tr; }


    Tr_vertex_handle const& center_vertex() const { return m_center_vertex; }
    Tr_vertex_handle & center_vertex() { return m_center_vertex; }

  private:
    Triangulation* m_tr;
    Tr_vertex_handle m_center_vertex;
  };

public:
  typedef typename std::vector<Tangent_space_basis>     TS_container;
  typedef typename std::vector<Orthogonal_space_basis>  OS_container;
private:
  typedef typename std::vector<Tr_and_VH>               Tr_container;
  typedef typename std::vector<Vector>                  Vectors;

  // An Incident_simplex is the list of the vertex indices
  // except the center vertex
  typedef std::set<std::size_t>                         Incident_simplex;
  typedef std::vector<Incident_simplex>                 Star;
  typedef std::vector<Star>                             Stars_container;

  // For transform_iterator
  static const Tr_point &vertex_handle_to_point(Tr_vertex_handle vh)
  {
    return vh->point();
  }
  template <typename P, typename VH>
  static const P &vertex_handle_to_point(VH vh)
  {
    return vh->point();
  }

public:
  typedef Tangential_complex_::Simplicial_complex     Simplicial_complex;

  /// Constructor for a range of points
  template <typename InputIterator>
  Tangential_complex(InputIterator first, InputIterator last,
                     double sparsity, int intrinsic_dimension,
#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
                     InputIterator first_for_tse, InputIterator last_for_tse,
#endif
                     const Kernel &k = Kernel()
                     )
  : m_k(k),
    m_intrinsic_dimension(intrinsic_dimension),
    m_half_sparsity(0.5*sparsity),
    m_sq_half_sparsity(m_half_sparsity*m_half_sparsity),
    m_ambient_dim(k.point_dimension_d_object()(*first)),
    m_points(first, last),
    m_weights(m_points.size(), FT(0))
#ifdef CGAL_TC_PERTURB_WEIGHT
    , m_weights_memory()
#endif
#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_PERTURB_POSITION) \
  && defined(CGAL_TC_GLOBAL_REFRESH)
    , m_p_perturb_mutexes(NULL)
#endif
    , m_points_ds(m_points)
    , m_are_tangent_spaces_computed(m_points.size(), false)
    , m_tangent_spaces(m_points.size(), Tangent_space_basis())
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
    , m_orth_spaces(m_points.size(), Orthogonal_space_basis())
#endif
#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
    , m_points_for_tse(first_for_tse, last_for_tse)
    , m_points_ds_for_tse(m_points_for_tse)
#endif
  { }

  /// Destructor
  ~Tangential_complex()
  {
#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_PERTURB_POSITION) \
 && defined(CGAL_TC_GLOBAL_REFRESH)
    delete [] m_p_perturb_mutexes;
#endif
  }

  int intrinsic_dimension() const
  {
    return m_intrinsic_dimension;
  }
  int ambient_dimension() const
  {
    return m_ambient_dim;
  }

  std::size_t number_of_vertices() const
  {
    return m_points.size();
  }

  void set_weights(const Weights& weights)
  {
    m_weights = weights;
#ifdef CGAL_TC_PERTURB_WEIGHT
    m_weights_memory = weights;
#endif
  }

  void set_tangent_planes(const TS_container& tangent_spaces
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
                        , const OS_container& orthogonal_spaces
#endif
                         )
  {
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    std::cerr << "Cannot use CGAL_TC_PERTURB_TANGENT_SPACE and set "
              << " tangent spaces manually at the same time" << std::endl;
    std::exit(EXIT_FAILURE);
#endif
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
    CGAL_assertion(m_points.size() == tangent_spaces.size()
                   && m_points.size() == orthogonal_spaces.size());
#else
    CGAL_assertion(m_points.size() == tangent_spaces.size());
#endif
    m_tangent_spaces = tangent_spaces;
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
    m_orth_spaces = orthogonal_spaces;
#endif
    for(std::size_t i=0; i<m_points.size(); ++i)
      m_are_tangent_spaces_computed[i] = true;
  }

  void compute_tangential_complex()
  {
#if defined(CGAL_TC_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    Wall_clock_timer t;
    ttt_intersect = 0;
    ttt_star = 0;
#endif

    // We need to do that because we don't want the container to copy the
    // already-computed triangulations (while resizing) since it would
    // invalidate the vertex handles stored beside the triangulations
    m_triangulations.resize(m_points.size());
    m_stars.resize(m_points.size());
#ifdef CGAL_LINKED_WITH_TBB
    //m_tr_mutexes.resize(m_points.size());
#endif
#ifdef CGAL_TC_PERTURB_POSITION
    m_translations.resize(m_points.size(),
                          m_k.construct_vector_d_object()(m_ambient_dim));
# if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
    delete [] m_p_perturb_mutexes;
    m_p_perturb_mutexes = new Mutex_for_perturb[m_points.size()];
# endif
#endif
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    m_perturb_tangent_space.resize(m_points.size(), false);
#endif

#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()),
        Compute_tangent_triangulation(*this)
      );
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      for (std::size_t i = 0 ; i < m_points.size() ; ++i)
#ifdef CGAL_ALPHA_TC
        compute_alpha_tangent_triangulation(i, ALPHA);
#else
        compute_tangent_triangulation(i);
#endif
    }

#if defined(CGAL_TC_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    std::cerr << "Tangential complex computed in " << t.elapsed()
              << " seconds." << std::endl;
    std::cerr << "Intersect: " << ((double)ttt_intersect)/1000000 << " s\n"
              << "Star: " << ((double)ttt_star)/1000000 << " s\n";
#endif
  }

  void estimate_intrinsic_dimension() const
  {
    // Kernel functors
    typename Kernel::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();

    std::vector<FT> sum_eigen_values(m_ambient_dim, FT(0));
    std::size_t num_points_for_pca =
      std::pow(BASE_VALUE_FOR_PCA, m_intrinsic_dimension);

    typename Points::const_iterator it_p = m_points.begin();
    typename Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for ( ; it_p != it_p_end ; ++it_p)
    {
      const Point &p = *it_p;

      KNS_range kns_range = m_points_ds.query_ANN(p, num_points_for_pca, false);
      //******************************* PCA *************************************

      // One row = one point
      Eigen::MatrixXd mat_points(num_points_for_pca, m_ambient_dim);
      KNS_iterator nn_it = kns_range.begin();
      for (int j = 0 ;
            j < num_points_for_pca && nn_it != kns_range.end() ;
            ++j, ++nn_it)
      {
        for (int i = 0 ; i < m_ambient_dim ; ++i)
          mat_points(j, i) = CGAL::to_double(coord(m_points[nn_it->first], i));
      }
      Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
      Eigen::MatrixXd cov = centered.adjoint() * centered;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

      // The eigenvectors are sorted in increasing order of their corresponding
      // eigenvalues
      for (int i = 0 ; i < m_ambient_dim ; ++i)
        sum_eigen_values[i] += eig.eigenvalues()[i];

      //*************************************************************************
    }

    // CJTODO: replace this by an actual estimation
    for (FT v : sum_eigen_values) // CJTODO C++11
    {
      std::cout << v << " ";
    }
    std::cout << "\n";
  }

  void refresh_tangential_complex()
  {
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif
#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()),
        Compute_tangent_triangulation(*this)
      );
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      for (std::size_t i = 0 ; i < m_points.size() ; ++i)
      {
#ifdef CGAL_ALPHA_TC
        compute_alpha_tangent_triangulation(i, ALPHA);
#else
        compute_tangent_triangulation(i);
#endif
      }
    }

#ifdef CGAL_TC_PROFILING
    std::cerr << "Tangential complex refreshed in " << t.elapsed()
              << " seconds." << std::endl;
#endif
  }

  // time_limit in seconds: 0 = no fix to do, < 0 = no time limit
  Fix_inconsistencies_status fix_inconsistencies_using_perturbation(
    unsigned int &num_steps,
    std::size_t &initial_num_inconsistent_local_tr,
    std::size_t &best_num_inconsistent_local_tr,
    std::size_t &final_num_inconsistent_local_tr,
    double time_limit = -1.)
  {
    if (time_limit == 0.)
      return TIME_LIMIT_REACHED;

    Wall_clock_timer t;

#ifdef CGAL_TC_VERBOSE
    std::cerr << "Fixing inconsistencies..." << std::endl;
#endif

#ifdef CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
    std::pair<std::size_t, std::size_t> stats_before =
      number_of_inconsistent_simplices(false);

# ifdef CGAL_TC_VERBOSE
      std::cerr << "Initial number of inconsistencies: "
      << stats_before.second << std::endl;
# endif

    if (stats_before.second == 0)
    {
# ifdef CGAL_TC_VERBOSE
      std::cerr << "Nothing to fix." << std::endl;
# endif
      return 0;
    }
#endif // CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

    bool done = false;
    best_num_inconsistent_local_tr = m_triangulations.size();
    num_steps = 0;
    while (!done)
    {
      std::size_t num_inconsistent_local_tr = 0;

#ifdef CGAL_TC_PROFILING
      Wall_clock_timer t_fix_step;
#endif

      // Parallel
#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
      if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
      {
        tbb::combinable<std::size_t> num_inconsistencies;
        tbb::parallel_for(
          tbb::blocked_range<size_t>(0, m_triangulations.size()),
          Try_to_solve_inconsistencies_in_a_local_triangulation(
            *this, num_inconsistencies)
        );
        num_inconsistent_local_tr =
          num_inconsistencies.combine(std::plus<std::size_t>());
      }
      // Sequential
      else
#endif // CGAL_LINKED_WITH_TBB
      {
        for (std::size_t i = 0 ; i < m_triangulations.size() ; ++i)
        {
          num_inconsistent_local_tr +=
            (try_to_solve_inconsistencies_in_a_local_triangulation(i) ? 1 : 0);
        }
      }

#ifdef CGAL_TC_PROFILING
      std::cerr << "Attempt to fix inconsistencies: " << t_fix_step.elapsed()
                << " seconds." << std::endl;
#endif

#ifdef CGAL_TC_GLOBAL_REFRESH
      refresh_tangential_complex();
#endif

#ifdef CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
      std::pair<std::size_t, std::size_t> stats_after =
        number_of_inconsistent_simplices(false);

      std::cerr << std::endl
        << "=========================================================="
        << std::endl
        << "Inconsistencies (detailed stats):\n"
        << "  * Number of vertices: " << m_points.size() << std::endl
        << std::endl
        << "  * BEFORE fix_inconsistencies_using_perturbation:" << std::endl
        << "    - Total number of simplices in stars (incl. duplicates): "
        << stats_before.first << std::endl
        << "    - Num inconsistent simplices in stars (incl. duplicates): "
        << stats_before.second
        << " (" << 100. * stats_before.second / stats_before.first << "%)"
        << std::endl
        << "  * Num inconsistent stars: "
        << num_inconsistent_local_tr
        << " (" << 100. * num_inconsistent_local_tr / m_points.size() << "%)"
        << std::endl
        << std::endl
        << "  * AFTER fix_inconsistencies_using_perturbation:" << std::endl
        << "    - Total number of simplices in stars (incl. duplicates): "
        << stats_after.first << std::endl
        << "    - Num inconsistent simplices in stars (incl. duplicates): "
        << stats_after.second << std::endl
        << "    - Percentage of inconsistencies: "
        << 100. * stats_after.second / stats_after.first << "%"
        << std::endl
        << "=========================================================="
        << std::endl;

      stats_before = stats_after;

#else // CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
# ifdef CGAL_TC_VERBOSE
      std::cerr << std::endl
        << "=========================================================="
        << std::endl
        << "fix_inconsistencies_using_perturbation():\n"
        << "  * " << m_points.size() << " vertices" << std::endl
        << "  * " << num_inconsistent_local_tr
        << " (" << 100. * num_inconsistent_local_tr / m_points.size() << "%)"
        << " inconsistent stars encountered" << std::endl
        << "=========================================================="
        << std::endl;
# endif
#endif // CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

      if (num_steps == 0)
        initial_num_inconsistent_local_tr = num_inconsistent_local_tr;

      if (num_inconsistent_local_tr < best_num_inconsistent_local_tr)
        best_num_inconsistent_local_tr = num_inconsistent_local_tr;

      final_num_inconsistent_local_tr = num_inconsistent_local_tr;

      ++num_steps;
      done = (num_inconsistent_local_tr == 0);
      if (!done && time_limit > 0. && t.elapsed() > time_limit)
      {
#ifdef CGAL_TC_VERBOSE
        std::cerr << "Time limit reached." << std::endl;
#endif
        return TIME_LIMIT_REACHED;
      }
    }

    return TC_FIXED;
  }

  // Return a pair<num_simplices, num_inconsistent_simplices>
  std::pair<std::size_t, std::size_t> number_of_inconsistent_simplices(
#ifdef CGAL_TC_VERBOSE
    bool verbose = true
#else
    bool verbose = false
#endif
    ) const
  {
    std::size_t num_simplices = 0;
    std::size_t num_inconsistent_simplices = 0;
    typename Tr_container::const_iterator it_tr = m_triangulations.begin();
    typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for (std::size_t idx = 0 ; it_tr != it_tr_end ; ++it_tr, ++idx)
    {
      // For each cell
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        // Don't export infinite cells
        if (*it_inc_simplex->rbegin() == std::numeric_limits<std::size_t>::max())
          continue;

        std::set<std::size_t> c = *it_inc_simplex;
        c.insert(idx); // Add the missing index

        if (!is_simplex_consistent(c))
          ++num_inconsistent_simplices;

        ++num_simplices;
      }
    }

    if (verbose)
    {
      std::cerr << std::endl
        << "=========================================================="
        << std::endl
        << "Inconsistencies:\n"
        << "  * Number of vertices: " << m_points.size() << std::endl
        << "  * Total number of simplices in stars (incl. duplicates): "
        << num_simplices << std::endl
        << "  * Number of inconsistent simplices in stars (incl. duplicates): "
        << num_inconsistent_simplices << std::endl
        << "  * Percentage of inconsistencies: "
        << 100 * num_inconsistent_simplices / num_simplices << "%" << std::endl
        << "=========================================================="
        << std::endl;
    }

    return std::make_pair(num_simplices, num_inconsistent_simplices);
  }

  // Return the max dimension of the simplices
  int export_TC(Simplicial_complex &complex,
    bool export_infinite_simplices = false) const
  {
    int max_dim = -1;

    typename Tr_container::const_iterator it_tr = m_triangulations.begin();
    typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for (std::size_t idx = 0 ; it_tr != it_tr_end ; ++it_tr, ++idx)
    {
      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        // Don't export infinite cells
        if (!export_infinite_simplices && *it_inc_simplex->rbegin()
            == std::numeric_limits<std::size_t>::max())
          continue;

        std::set<std::size_t> c = *it_inc_simplex;
        if (static_cast<int>(c.size()) > max_dim)
          max_dim = static_cast<int>(c.size());
        // Add the missing center vertex
        c.insert(idx);
        complex.add_simplex(c);
      }
    }
    return max_dim;
  }

  void check_and_solve_inconsistencies_by_adding_higher_dim_simplices()
  {
    // CJTODO: parallel_for???
    for (std::size_t idx = 0 ; idx < m_triangulations.size() ; ++idx)
    {
      bool inconsistencies_found = false;
      do
      {
        Star::const_iterator it_inc_simplex = m_stars[idx].begin();
        Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
        for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
        {
          inconsistencies_found =
            check_and_solve_inconsistencies_by_adding_higher_dim_simplices(
              idx, *it_inc_simplex);

          // m_stars[idx] has been modified, let's start again
          // CJTODO: optimize?
          if (inconsistencies_found)
            break;
        }
      } while (inconsistencies_found);
    }

    // CJTODO TEMP
    std::pair<std::size_t, std::size_t> stats_after =
      number_of_inconsistent_simplices(false);
    std::cerr << "AFTER check_and_solve_inconsistencies_by_adding_higher_dim_simplices():\n"
      << "    - Total number of simplices in stars (incl. duplicates): "
      << stats_after.first << std::endl
      << "    - Num inconsistent simplices in stars (incl. duplicates): "
      << stats_after.second << std::endl
      << "    - Percentage of inconsistencies: "
      << 100. * stats_after.second / stats_after.first << "%"
      << std::endl;
  }

  std::ostream &export_to_off(
    const Simplicial_complex &complex, std::ostream & os,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_red = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_green = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_blue = NULL)
    const
  {
    return export_to_off(
      os, false, p_simpl_to_color_in_red, p_simpl_to_color_in_green, 
      p_simpl_to_color_in_blue, &complex);
  }

  std::ostream &export_to_off(
    std::ostream & os, bool color_inconsistencies = false,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_red = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_green = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_blue = NULL,
    const Simplicial_complex *p_complex = NULL) const
  {
    if (m_points.empty())
      return os;

    if (m_ambient_dim < 2)
    {
      std::cerr << "Error: export_to_off => ambient dimension should be >= 2."
                << std::endl;
      os << "Error: export_to_off => ambient dimension should be >= 2."
         << std::endl;
      return os;
    }
    if (m_ambient_dim > 3)
    {
      std::cerr << "Warning: export_to_off => ambient dimension should be "
                   "<= 3. Only the first 3 coordinates will be exported."
                << std::endl;
    }

    if (m_intrinsic_dimension < 1 || m_intrinsic_dimension > 3)
    {
      std::cerr << "Error: export_to_off => intrinsic dimension should be "
                   "between 1 and 3."
                << std::endl;
      os << "Error: export_to_off => intrinsic dimension should be "
            "between 1 and 3."
         << std::endl;
      return os;
    }

    std::stringstream output;
    std::size_t num_simplices, num_vertices;
    export_vertices_to_off(output, num_vertices);
    if (p_complex)
    {
      export_simplices_to_off(
        *p_complex, output, num_simplices, p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, p_simpl_to_color_in_blue);
    }
    else
    {
      export_simplices_to_off(
        output, num_simplices, color_inconsistencies, p_simpl_to_color_in_red,
        p_simpl_to_color_in_green, p_simpl_to_color_in_blue);
    }

#ifdef CGAL_TC_EXPORT_NORMALS
    os << "N";
#endif

    os << "OFF \n"
       << num_vertices << " "
       << num_simplices << " "
       << "0 \n"
       << output.str();

    return os;
  }
  
  // Return a pair<num_simplices, num_inconsistent_simplices>
  void export_inconsistent_stars_to_OFF_files(
    std::string const& filename_base) const
  {
    std::size_t num_simplices = 0;
    std::size_t num_inconsistent_simplices = 0;
    typename Tr_container::const_iterator it_tr = m_triangulations.begin();
    typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for (std::size_t idx = 0 ; it_tr != it_tr_end ; ++it_tr, ++idx)
    {
      // We build a SC along the way in case it's inconsistent
      Simplicial_complex sc;
      // For each cell
      bool is_inconsistent = false;
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; 
           ++it_inc_simplex)
      {
        // Skip infinite cells
        if (*it_inc_simplex->rbegin() == std::numeric_limits<std::size_t>::max())
          continue;

        std::set<std::size_t> c = *it_inc_simplex;
        c.insert(idx); // Add the missing index

        sc.add_simplex(c);

        // If we do not already know this star is inconsistent, test it
        if (!is_inconsistent && !is_simplex_consistent(c))
          is_inconsistent = true;
      }

      if (is_inconsistent)
      {
        // Export star to OFF file
        std::stringstream output_filename;
        output_filename << filename_base << "_" << idx << ".off";
        std::ofstream off_stream(output_filename.str().c_str());
        export_to_off(sc, off_stream);
      }
    }
  }
  

  bool check_if_all_simplices_are_in_the_ambient_delaunay(
    const Simplicial_complex *p_complex = NULL,
    bool check_for_any_dimension_simplices = true,
    std::set<std::set<std::size_t> > * incorrect_simplices = NULL) const
  {
    typedef Simplicial_complex::Simplex                     Simplex;
    typedef Simplicial_complex::Simplex_range               Simplex_range;

    if (m_points.empty())
      return true;

    typedef Regular_triangulation_euclidean_traits<Kernel>    RT_Traits;
    typedef Regular_triangulation<
      RT_Traits,
      Triangulation_data_structure<
        typename RT_Traits::Dimension,
        Triangulation_vertex<RT_Traits, Vertex_data>
      > >                                                     RT;
    typedef typename RT::Vertex_handle                        RT_VH;
    typedef typename RT::Finite_full_cell_const_iterator      FFC_it;

    //-------------------------------------------------------------------------
    // Build the ambient Delaunay triangulation
    // Then save its simplices into "amb_dt_simplices"
    //-------------------------------------------------------------------------

    RT ambient_dt(m_ambient_dim);
    for (std::size_t i=0; i<m_points.size(); ++i)
    {
      const Weighted_point wp = compute_perturbed_weighted_point(i);
      RT_VH vh = ambient_dt.insert(wp);
      vh->data() = i;
    }

    std::set<Simplex> amb_dt_simplices;

    for (FFC_it cit = ambient_dt.finite_full_cells_begin() ;
         cit != ambient_dt.finite_full_cells_end() ; ++cit )
    {
      int lowest_dim =
        (check_for_any_dimension_simplices ? 1 : m_intrinsic_dimension);
      int highest_dim =
        (check_for_any_dimension_simplices ? m_ambient_dim : m_intrinsic_dimension);

      for (int dim = lowest_dim ; dim <= highest_dim ; ++dim)
      {
        CGAL::Combination_enumerator<int> combi(dim + 1, 0, m_ambient_dim + 1);

        for ( ; !combi.finished() ; ++combi)
        {
          Simplex simplex;
          for (int i = 0 ; i < dim + 1 ; ++i)
            simplex.insert(cit->vertex(combi[i])->data());

          amb_dt_simplices.insert(simplex);
        }
      }
    }

    //-------------------------------------------------------------------------
    // If p_complex is NULL, parse the TC and
    // save its simplices into "stars_simplices"
    //-------------------------------------------------------------------------

    Simplex_range const *p_simplices;

    if (!p_complex)
    {
      Simplex_range stars_simplices;

      typename Tr_container::const_iterator it_tr = m_triangulations.begin();
      typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
      // For each triangulation
      for ( ; it_tr != it_tr_end ; ++it_tr)
      {
        Triangulation const& tr    = it_tr->tr();
        Tr_vertex_handle center_vh = it_tr->center_vertex();

        std::vector<Tr_full_cell_handle> incident_cells;
        tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

        typename std::vector<Tr_full_cell_handle>::const_iterator it_c =
                                                           incident_cells.begin();
        typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end =
                                                             incident_cells.end();
        // For each cell
        for ( ; it_c != it_c_end ; ++it_c)
        {
          if (tr.is_infinite(*it_c))
          {
            std::cerr << "Warning: infinite cell in star" << std::endl;
            continue;
          }
          Simplex simplex;
          for (int i = 0 ; i < tr.current_dimension() + 1 ; ++i)
            simplex.insert((*it_c)->vertex(i)->data());

          stars_simplices.insert(simplex);
        }
      }

      p_simplices = &stars_simplices;
    }
    else
    {
      p_simplices = &p_complex->simplex_range();
    }

    //-------------------------------------------------------------------------
    // Check if simplices of "*p_complex" are all in "amb_dt_simplices"
    //-------------------------------------------------------------------------

    std::set<Simplex> diff;
    if (!incorrect_simplices)
      incorrect_simplices = &diff;
    std::set_difference(p_simplices->begin(),  p_simplices->end(),
                   amb_dt_simplices.begin(), amb_dt_simplices.end(),
                   std::inserter(*incorrect_simplices,
                                 incorrect_simplices->begin()) );

    if (!incorrect_simplices->empty())
    {
      std::cerr
        << "ERROR check_if_all_simplices_are_in_the_ambient_delaunay:"
        << std::endl
        << "  Number of simplices in ambient RT: " << amb_dt_simplices.size()
        << std::endl
        << "  Number of unique simplices in TC stars: " << p_simplices->size()
        << std::endl
        << "  Number of wrong simplices: " << incorrect_simplices->size()
        << std::endl;
      return false;
    }
    else
    {
#ifdef CGAL_TC_VERBOSE
      std::cerr
        << "SUCCESS check_if_all_simplices_are_in_the_ambient_delaunay:"
        << std::endl
        << "  Number of simplices in ambient RT: " << amb_dt_simplices.size()
        << std::endl
        << "  Number of unique simplices in TC stars: " << p_simplices->size()
        << std::endl
        << "  Number of wrong simplices: " << incorrect_simplices->size()
        << std::endl;
#endif
      return true;
    }
  }

private:

  class Compare_distance_to_ref_point
  {
  public:
    Compare_distance_to_ref_point(Point const& ref, Kernel const& k)
      : m_ref(ref), m_k(k) {}

    bool operator()(Point const& p1, Point const& p2)
    {
      typename Kernel::Squared_distance_d sqdist =
        m_k.squared_distance_d_object();
      return sqdist(p1, m_ref) < sqdist(p2, m_ref);
    }

  private:
    Point const& m_ref;
    Kernel const& m_k;
  };

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for compute_tangential_complex function
  class Compute_tangent_triangulation
  {
    Tangential_complex & m_tc;

  public:
    // Constructor
    Compute_tangent_triangulation(
      Tangential_complex &tc)
    : m_tc(tc)
    { }

    // Constructor
    Compute_tangent_triangulation(const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc)
    { }

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
      {
#ifdef CGAL_ALPHA_TC
        m_tc.compute_alpha_tangent_triangulation(i, ALPHA);
#else
        m_tc.compute_tangent_triangulation(i);
#endif
      }
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  void compute_tangent_triangulation(std::size_t i, bool verbose = false)
  {
    if (verbose)
      std::cerr << "** Computing tangent tri #" << i << " **" << std::endl;
    //std::cerr << "***********************************************" << std::endl;
    Triangulation &local_tr =
      m_triangulations[i].construct_triangulation(m_intrinsic_dimension);
    const Tr_traits &local_tr_traits = local_tr.geom_traits();
    Tr_vertex_handle &center_vertex = m_triangulations[i].center_vertex();

    // Kernel functor & objects
    typename Kernel::Squared_distance_d k_sqdist =
      m_k.squared_distance_d_object();

    // Triangulation's traits functor & objects
    typename Tr_traits::Point_weight_d point_weight =
      local_tr_traits.point_weight_d_object();
    typename Tr_traits::Power_center_d power_center =
      local_tr_traits.power_center_d_object();

    // No need to lock the mutex here since this will not be called while
    // other threads are perturbing the positions
    const Point center_pt = compute_perturbed_point(i);

    // Estimate the tangent space
    if (!m_are_tangent_spaces_computed[i])
    {
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
      m_tangent_spaces[i] =
        compute_tangent_space(center_pt, i, true/*normalize*/, &m_orth_spaces[i]);
#else
      m_tangent_spaces[i] = compute_tangent_space(center_pt, i);
#endif
    }
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    else if (m_perturb_tangent_space[i])
    {
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
      m_tangent_spaces[i] = compute_tangent_space(center_pt, i,
                                                  true /*normalize_basis*/,
                                                  &m_orth_spaces[i],
                                                  true /*perturb*/);
#else
      m_tangent_spaces[i] = compute_tangent_space(center_pt, i,
                                                  true /*normalize_basis*/,
                                                  NULL /*ortho basis*/,
                                                  true /*perturb*/);
#endif
      m_perturb_tangent_space[i] = false;
    }
#endif

    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p)
    //***************************************************

    // Insert p
    typename Kernel::Equal_d eq = m_k.equal_d_object();

    Tr_point proj_wp;
    if(eq(compute_perturbed_point(i), m_tangent_spaces[i].origin()))
    {
      proj_wp = local_tr_traits.construct_weighted_point_d_object()(
      local_tr_traits.construct_point_d_object()(m_intrinsic_dimension, ORIGIN),
      m_weights[i]);
    }
    else
    {
      const Weighted_point& wp = compute_perturbed_weighted_point(i);
      proj_wp = project_point_and_compute_weight(wp, m_tangent_spaces[i],
                                                 local_tr_traits);
    }

    center_vertex = local_tr.insert(proj_wp);
    center_vertex->data() = i;
    if (verbose)
      std::cerr << "* Inserted point #" << i << std::endl;

    //const int NUM_NEIGHBORS = 150;
    //KNS_range ins_range = m_points_ds.query_ANN(center_pt, NUM_NEIGHBORS);
    INS_range ins_range = m_points_ds.query_incremental_ANN(center_pt);

    // While building the local triangulation, we keep the radius
    // of the sphere "star sphere" centered at "center_vertex"
    // and which contains all the
    // circumspheres of the star of "center_vertex"
    boost::optional<FT> squared_star_sphere_radius_plus_margin;

    // Insert points until we find a point which is outside "star shere"
    for (INS_iterator nn_it = ins_range.begin() ;
         nn_it != ins_range.end() ;
         ++nn_it)
    {
      std::size_t neighbor_point_idx = nn_it->first;

      // ith point = p, which is already inserted
      if (neighbor_point_idx != i)
      {
        // No need to lock the Mutex_for_perturb here since this will not be
        // called while other threads are perturbing the positions
        Point neighbor_pt;
        FT neighbor_weight;
        compute_perturbed_weighted_point(
          neighbor_point_idx, neighbor_pt, neighbor_weight);

        // "4*m_sq_half_sparsity" because both points can be perturbed
        if (squared_star_sphere_radius_plus_margin
          && k_sqdist(center_pt, neighbor_pt)
             > *squared_star_sphere_radius_plus_margin)
          break;

        Tr_point proj_pt = project_point_and_compute_weight(
          neighbor_pt, neighbor_weight, m_tangent_spaces[i],
          local_tr_traits);

        Tr_vertex_handle vh = local_tr.insert_if_in_star(proj_pt, center_vertex);
        //Tr_vertex_handle vh = local_tr.insert(proj_pt);
        if (vh != Tr_vertex_handle())
        {
          if (verbose)
            std::cerr << "* Inserted point #" << neighbor_point_idx << std::endl;

          vh->data() = neighbor_point_idx;

          // Let's recompute squared_star_sphere_radius_plus_margin
          if (local_tr.current_dimension() >= m_intrinsic_dimension)
          {
            squared_star_sphere_radius_plus_margin = boost::none;
            // Get the incident cells and look for the biggest circumsphere
            std::vector<Tr_full_cell_handle> incident_cells;
            local_tr.incident_full_cells(
              center_vertex,
              std::back_inserter(incident_cells));
            for (typename std::vector<Tr_full_cell_handle>::iterator cit =
                 incident_cells.begin(); cit != incident_cells.end(); ++cit)
            {
              Tr_full_cell_handle cell = *cit;
              if (local_tr.is_infinite(cell))
              {
                squared_star_sphere_radius_plus_margin = boost::none;
                break;
              }
              else
              {
                Tr_point c = power_center(
                  boost::make_transform_iterator(
                    cell->vertices_begin(),
                    vertex_handle_to_point<Tr_point, Tr_vertex_handle>),
                  boost::make_transform_iterator(
                    cell->vertices_end(),
                    vertex_handle_to_point<Tr_point, Tr_vertex_handle>));

                FT sq_power_sphere_diam = 4*point_weight(c);

                if (!squared_star_sphere_radius_plus_margin
                 || sq_power_sphere_diam >
                    *squared_star_sphere_radius_plus_margin)
                {
                  squared_star_sphere_radius_plus_margin = sq_power_sphere_diam;
                }
              }
            }

            // Let's add the margin, now
            // The value depends on whether we perturb weight or position
            if (squared_star_sphere_radius_plus_margin)
            {
#ifdef CGAL_TC_PERTURB_WEIGHT
              squared_star_sphere_radius_plus_margin =
                *squared_star_sphere_radius_plus_margin + 4*m_sq_half_sparsity;
#else
              squared_star_sphere_radius_plus_margin = CGAL::square(
                CGAL::sqrt(*squared_star_sphere_radius_plus_margin)
                + 2*m_half_sparsity);
#endif
            }
          }
        }
      }
    }

    //***************************************************
    // Update the associated star (in m_stars)
    //***************************************************
    Star &star = m_stars[i];
    star.clear();
    int cur_dim_plus_1 = local_tr.current_dimension() + 1;

    std::vector<Tr_full_cell_handle> incident_cells;
    local_tr.incident_full_cells(
      center_vertex, std::back_inserter(incident_cells));

    typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
    typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end= incident_cells.end();
    // For each cell
    for ( ; it_c != it_c_end ; ++it_c)
    {
      // Will contain all indices except center_vertex
      Incident_simplex incident_simplex;
      for (int j = 0 ; j < cur_dim_plus_1 ; ++j)
      {
        std::size_t index = (*it_c)->vertex(j)->data();
        if (index != i)
          incident_simplex.insert(index);
      }
      star.push_back(incident_simplex);
    }

    // CJTODO DEBUG
    //std::cerr << "\nChecking topology and geometry..."
    //          << (local_tr.is_valid(true) ? "OK.\n" : "Error.\n");
    // DEBUG: output the local mesh into an OFF file
    //std::stringstream sstr;
    //sstr << "data/local_tri_" << i << ".off";
    //std::ofstream off_stream_tr(sstr.str());
    //CGAL::export_triangulation_to_off(off_stream_tr, local_tr);
  }

#ifdef CGAL_ALPHA_TC
  void compute_alpha_tangent_triangulation(std::size_t i, FT alpha,
                                           bool verbose = false)
  {
    if (verbose)
      std::cerr << "** Computing alpha tangent tri #"
                << i << " **" << std::endl;

    typedef Regular_triangulation_euclidean_traits<Kernel>    Amb_RT_Traits;
    typedef Regular_triangulation<
      Amb_RT_Traits,
      Triangulation_data_structure<
        typename Amb_RT_Traits::Dimension,
        Triangulation_vertex<Amb_RT_Traits, Vertex_data>
      > >                                                     Amb_RT;
    typedef typename Amb_RT::Weighted_point                   Amb_RT_Point;
    typedef typename Amb_RT::Vertex_handle                    Amb_RT_VH;
    typedef typename Amb_RT::Full_cell_handle                 Amb_RT_FCH;
    //typedef typename Amb_RT::Finite_full_cell_const_iterator  Amb_FFC_it;

    typename Kernel::Point_drop_weight_d k_drop_w =
      m_k.point_drop_weight_d_object();
    typename Kernel::Squared_distance_d k_sqdist =
      m_k.squared_distance_d_object();
    typename Kernel::Point_weight_d k_point_weight =
      m_k.point_weight_d_object();
    typename Kernel::Power_center_d k_power_center =
      m_k.power_center_d_object();

    // No need to lock the mutex here since this will not be called while
    // other threads are perturbing the positions
    const Point &center_pt = m_points[i];

    // Estimate the tangent space
    if (!m_are_tangent_spaces_computed[i])
    {
      m_tangent_spaces[i] =
        compute_tangent_space(center_pt, i, true /*normalize*/, &m_orth_spaces[i]);
    }
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    else if (m_perturb_tangent_space[i])
    {
      m_tangent_spaces[i] =
        compute_tangent_space(center_pt, i, true, &m_orth_spaces[i], true);
      m_perturb_tangent_space[i] = false;
    }
#endif

    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p in the AMBIENT triangulation)
    //***************************************************

#ifdef CGAL_TC_PROFILING
        Wall_clock_timer t_star;
#endif

    Amb_RT local_amb_tr(m_ambient_dim);

    // Insert p
    Weighted_point wp = m_k.construct_weighted_point_d_object()(center_pt,
                                                                m_weights[i]);
    Amb_RT_VH center_vertex = local_amb_tr.insert(wp);
    center_vertex->data() = i;
    if (verbose)
      std::cerr << "* Inserted point #" << i << std::endl;

    INS_range ins_range = m_points_ds.query_incremental_ANN(center_pt);

    // While building the local triangulation, we keep the radius
    // of the sphere "star sphere" centered at "center_vertex"
    // and which contains all the
    // circumspheres of the star of "center_vertex"
    boost::optional<FT> squared_star_sphere_radius_plus_margin;

    // Insert points until we find a point which is outside "star shere"
    for (INS_iterator nn_it = ins_range.begin() ;
         nn_it != ins_range.end() ;
         ++nn_it)
    {
      std::size_t neighbor_point_idx = nn_it->first;

      // ith point = p, which is already inserted
      if (neighbor_point_idx != i)
      {
        // No need to lock the Mutex_for_perturb here since this will not be
        // called while other threads are perturbing the positions
        Weighted_point neighbor_wp =
          compute_perturbed_weighted_point(neighbor_point_idx);
// fixme ?
// Above seems incorrect since we later pass m_points[] to
// "does_voronoi_face_and_alpha_tangent_subspace_intersect()" and we lose
// the possible position pertubation of wp.
// Either we ignore the translations and it's neighbor_wp = points[idx]+weight
// or need to give m_translations to the voronoi intersections computations

        // "4*m_sq_half_sparsity" because both points can be perturbed
        if (squared_star_sphere_radius_plus_margin
          && k_sqdist(center_pt, k_drop_w(neighbor_wp))
             > *squared_star_sphere_radius_plus_margin)
          break;

        Amb_RT_VH vh = local_amb_tr.insert_if_in_star(neighbor_wp, center_vertex);
        //Amb_RT_VH vh = local_amb_tr.insert(neighbor_wp);
        if (vh != Amb_RT_VH())
        {
          if (verbose)
            std::cerr << "* Inserted point #" << neighbor_point_idx << std::endl;

          vh->data() = neighbor_point_idx;

          // Let's recompute squared_star_sphere_radius_plus_margin
          if (local_amb_tr.current_dimension() == m_ambient_dim)
          {
            squared_star_sphere_radius_plus_margin = boost::none;
            // Get the incident cells and look for the biggest circumsphere
            std::vector<Amb_RT_FCH> incident_cells;
            local_amb_tr.incident_full_cells(
              center_vertex,
              std::back_inserter(incident_cells));
            for (typename std::vector<Amb_RT_FCH>::iterator cit =
                 incident_cells.begin(); cit != incident_cells.end(); ++cit)
            {
              Amb_RT_FCH cell = *cit;
              if (local_amb_tr.is_infinite(cell))
              {
                squared_star_sphere_radius_plus_margin = boost::none;
                break;
              }
              else
              {
                Amb_RT_Point c = k_power_center(
                  boost::make_transform_iterator(
                    cell->vertices_begin(),
                    vertex_handle_to_point<Weighted_point, Amb_RT_VH>),
                  boost::make_transform_iterator(
                    cell->vertices_end(),
                    vertex_handle_to_point<Weighted_point, Amb_RT_VH>));

                FT sq_power_sphere_diam = 4*k_point_weight(c);

                if (!squared_star_sphere_radius_plus_margin
                 || sq_power_sphere_diam >
                    *squared_star_sphere_radius_plus_margin)
                {
                  squared_star_sphere_radius_plus_margin = sq_power_sphere_diam;
                }
              }
            }

            // Let's add the margin, now
            // The value depends on whether we perturb weight or position
            if (squared_star_sphere_radius_plus_margin)
            {
#ifdef CGAL_TC_PERTURB_WEIGHT
              squared_star_sphere_radius_plus_margin =
                *squared_star_sphere_radius_plus_margin + 4*m_sq_half_sparsity;
#else
              squared_star_sphere_radius_plus_margin = CGAL::square(
                CGAL::sqrt(*squared_star_sphere_radius_plus_margin)
                + 2*m_half_sparsity);
#endif
            }
          }
        }
      }
    }

#if defined(CGAL_TC_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    ttt_star += 1000000*t_star.elapsed();
#endif

    //***************************************************
    // Parse the faces of the star and add the ones that are in the
    // restriction to alpha-Tp
    // Update the associated star (in m_stars)
    //***************************************************
    Star &star = m_stars[i];
    star.clear();
    int cur_dim_plus_1 = m_ambient_dim + 1;

    std::vector<Amb_RT_FCH> incident_cells;
    local_amb_tr.incident_full_cells(
      center_vertex, std::back_inserter(incident_cells));

    typedef std::set<std::size_t> DT_face; // DT face without center vertex (i)
    typedef std::set<std::size_t> Neighbor_vertices;
    typedef std::map<DT_face, Neighbor_vertices> DT_faces_and_neighbors;

    // Maps that associate a k-face and the list of its neighbor points
    // (i.e. there are k+1-cofaces that contain these points)
    // N.B.: each k-face contains 'i', so 'i' is not stored in the faces
    // faces_and_neighbors[0] => dim 1, faces_and_neighbors[1] => dim 2
    std::vector<DT_faces_and_neighbors> faces_and_neighbors;
    faces_and_neighbors.resize(m_ambient_dim);

    // Fill faces_and_neighbors
    // Let's first take care of the D-faces
    typename std::vector<Amb_RT_FCH>::const_iterator it_c = incident_cells.begin();
    typename std::vector<Amb_RT_FCH>::const_iterator it_c_end = incident_cells.end();
    // For each cell
    for ( ; it_c != it_c_end ; ++it_c)
    {
      DT_face face;
      // CJTODO: use (*it_c)->vertices_begin(), etc.
      for (int j = 0 ; j < cur_dim_plus_1 ; ++j)
      {
        std::size_t index = (*it_c)->vertex(j)->data();
        if (index == std::numeric_limits<std::size_t>::max())
          goto next_face;
        if (index != i)
          face.insert(index);
      }
      faces_and_neighbors[m_ambient_dim-1][face] = Neighbor_vertices();
next_face:
      ;
    }
    // Then the D-k-faces...
    int current_dim = m_ambient_dim - 1;
    while (current_dim > 0)
    {
      // Let's fill faces_and_neighbors[current_dim-1]
      // (stores the current_dim-faces)
      DT_faces_and_neighbors& cur_faces_and_nghb =
        faces_and_neighbors[current_dim-1];

      typedef DT_faces_and_neighbors::const_iterator FaN_it;
      // Parse k+1-faces
      for (FaN_it it_k_p1_face = faces_and_neighbors[current_dim].begin(),
                  it_k_p1_face_end = faces_and_neighbors[current_dim].end() ;
           it_k_p1_face != it_k_p1_face_end ; ++it_k_p1_face)
      {
        DT_face const& k_p1_face = it_k_p1_face->first;

        // Add each k faces to cur_faces_and_nghb
        std::size_t n = current_dim + 1; // Not +2 since 'i' is not stored
        std::vector<bool> booleans(n, false);
        std::fill(booleans.begin() + 1, booleans.end(), true);
        do
        {
          DT_face k_face;
          std::size_t remaining_vertex;
          DT_face::const_iterator it_v = k_p1_face.begin();
          for (std::size_t i = 0 ; i < n ; ++i, ++it_v)
          {
            if (booleans[i])
              k_face.insert(*it_v);
            else
              remaining_vertex = *it_v;
          }

          cur_faces_and_nghb[k_face].insert(remaining_vertex);

        } while (std::next_permutation(booleans.begin(), booleans.end()));
      }
      --current_dim;
    }

    // For each face V of Voronoi_cell(P[i]) - dim 0 to dim D-1
    // I.e. For each DT face F of the star - dim D to dim 1
    current_dim = m_ambient_dim;
    while (current_dim > 0)
    {
      // Remember: faces_and_neighbors[current_dim-1] stores
      // the current_dim-faces
      DT_faces_and_neighbors const& cur_faces_and_nghb =
        faces_and_neighbors[current_dim-1];

      for (DT_faces_and_neighbors::const_iterator
             it_f = cur_faces_and_nghb.begin(),
             it_f_end = cur_faces_and_nghb.end() ;
           it_f != it_f_end ; ++it_f)
      {
        Neighbor_vertices const& curr_neighbors = it_f->second;

        DT_face const& current_DT_face = it_f->first;
        CGAL_assertion(static_cast<int>(current_DT_face.size())
                       == current_dim);

        // P: list of current_DT_face points (including 'i')
        std::vector<std::size_t> P(
          current_DT_face.begin(), current_DT_face.end());
        P.push_back(i);

#ifdef CGAL_TC_PROFILING
        Wall_clock_timer t_inters;
#endif
        bool does_intersect =
          does_voronoi_face_and_alpha_tangent_subspace_intersect(
              m_points, m_weights, i, P, curr_neighbors,
              m_orth_spaces[i], alpha, m_k);
#if defined(CGAL_TC_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
        ttt_intersect += 1000000*t_inters.elapsed();
#endif
        if (does_intersect)
        {
          star.push_back(current_DT_face);

          // Clear all subfaces of current_DT_face from the maps
          for (int dim = current_dim - 1 ; dim > 0 ; --dim)
          {
            std::size_t n = current_DT_face.size();
            std::vector<bool> booleans(n, false);
            std::fill(booleans.begin() + n - dim, booleans.end(), true);
            do
            {
              DT_face dim_face;
              DT_face::const_iterator it_v = current_DT_face.begin();
              for (std::size_t i = 0 ; i < n ; ++i, ++it_v)
              {
                if (booleans[i])
                  dim_face.insert(*it_v);
              }

              faces_and_neighbors[dim-1].erase(dim_face);

            } while (std::next_permutation(booleans.begin(), booleans.end()));
          }
        }
      }

      --current_dim;
    }

    // CJTODO DEBUG
    //std::cerr << "\nChecking topology and geometry..."
    //          << (local_amb_tr.is_valid(true) ? "OK.\n" : "Error.\n");
    // DEBUG: output the local mesh into an OFF file
    //std::stringstream sstr;
    //sstr << "data/local_tri_" << i << ".off";
    //std::ofstream off_stream_tr(sstr.str());
    //CGAL::export_triangulation_to_off(off_stream_tr, local_amb_tr);
  }
#endif // CGAL_ALPHA_TC

  Tangent_space_basis compute_tangent_space(
      const Point &p
    , const std::size_t i
    , bool normalize_basis = true
    , Orthogonal_space_basis *p_orth_space_basis = NULL
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    , bool perturb = false
#endif
    )
  {
#ifdef CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_SPHERE_3

    // CJTODO: this is only for a sphere in R^3
    double tt1[3] = {-p[1] - p[2], p[0], p[0]};
    double tt2[3] = {p[1] * tt1[2] - p[2] * tt1[1],
                     p[2] * tt1[0] - p[0] * tt1[2],
                     p[0] * tt1[1] - p[1] * tt1[0]};
    Vector t1(3, &tt1[0], &tt1[3]);
    Vector t2(3, &tt2[0], &tt2[3]);

    // Normalize t1 and t2
    typename Kernel::Squared_length_d sqlen      = m_k.squared_length_d_object();
    typename Kernel::Scaled_vector_d  scaled_vec = m_k.scaled_vector_d_object();

    Tangent_space_basis ts;
    ts.reserve(m_intrinsic_dimension);
    ts.push_back(scaled_vec(t1, FT(1)/CGAL::sqrt(sqlen(t1))));
    ts.push_back(scaled_vec(t2, FT(1)/CGAL::sqrt(sqlen(t2))));

    m_are_tangent_spaces_computed[i] = true;

    return ts;

#elif defined(CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_TORUS_D)

    // CJTODO: this is only for torus_d
    Tangent_space_basis ts(p);
    ts.reserve(m_intrinsic_dimension);
    for (int dim = 0 ; dim < m_intrinsic_dimension ; ++dim)
    {
      std::vector<FT> tt(m_ambient_dim, 0.);
      tt[2*dim] = -p[2*dim + 1];
      tt[2*dim + 1] = p[2*dim];
      Vector t(2*m_intrinsic_dimension, tt.begin(), tt.end());
      ts.push_back(t);
    }

    m_are_tangent_spaces_computed[i] = true;
    
    //return compute_gram_schmidt_basis(ts, m_k);
    return ts;
    //******************************* PCA *************************************
    
#else

    unsigned int num_points_for_pca = static_cast<unsigned int>(
      std::pow(BASE_VALUE_FOR_PCA, m_intrinsic_dimension));

    // Kernel functors
    typename Kernel::Construct_vector_d      constr_vec =
      m_k.construct_vector_d_object();
    typename Kernel::Compute_coordinate_d    coord =
      m_k.compute_coordinate_d_object();
    typename Kernel::Squared_length_d        sqlen =
      m_k.squared_length_d_object();
    typename Kernel::Scaled_vector_d         scaled_vec =
      m_k.scaled_vector_d_object();
    typename Kernel::Scalar_product_d        inner_pdct =
      m_k.scalar_product_d_object();
    typename Kernel::Difference_of_vectors_d diff_vec =
      m_k.difference_of_vectors_d_object();
    //typename Kernel::Translated_point_d      transl =
    //  m_k.translated_point_d_object();

#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
    KNS_range kns_range = m_points_ds_for_tse.query_ANN(
      p, num_points_for_pca, false);
    const Points &points_for_pca = m_points_for_tse;
#else
    KNS_range kns_range = m_points_ds.query_ANN(p, num_points_for_pca, false);
    const Points &points_for_pca = m_points;
#endif

    // One row = one point
    Eigen::MatrixXd mat_points(num_points_for_pca, m_ambient_dim);
    KNS_iterator nn_it = kns_range.begin();
    for (unsigned int j = 0 ;
         j < num_points_for_pca && nn_it != kns_range.end() ;
         ++j, ++nn_it)
    {
      for (int i = 0 ; i < m_ambient_dim ; ++i)
      {
        //const Point p = transl(
        //  m_points[nn_it->first], m_translations[nn_it->first]);
        mat_points(j, i) = CGAL::to_double(coord(m_points[nn_it->first], i));
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
        if (perturb)
          mat_points(j, i) += m_random_generator.get_double(
            -m_half_sparsity, m_half_sparsity);
#endif
      }
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    Tangent_space_basis tsb(p); // p = compute_perturbed_point(i) here

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    for (int j = m_ambient_dim - 1 ;
         j >= m_ambient_dim - m_intrinsic_dimension ;
         --j)
    {
      if (normalize_basis)
      {
        Vector v = constr_vec(m_ambient_dim,
                              eig.eigenvectors().col(j).data(),
                              eig.eigenvectors().col(j).data() + m_ambient_dim);
        tsb.push_back(normalize_vector(v, m_k));
      }
      else
      {
        tsb.push_back(constr_vec(
          m_ambient_dim,
          eig.eigenvectors().col(j).data(),
          eig.eigenvectors().col(j).data() + m_ambient_dim));
      }
    }

    if (p_orth_space_basis)
    {
      p_orth_space_basis->origin() = p;
      for (int j = m_ambient_dim - m_intrinsic_dimension - 1 ;
           j >= 0 ;
           --j)
      {
        if (normalize_basis)
        {
          Vector v = constr_vec(m_ambient_dim,
                                eig.eigenvectors().col(j).data(),
                                eig.eigenvectors().col(j).data() + m_ambient_dim);
          p_orth_space_basis->push_back(normalize_vector(v, m_k));
        }
        else
        {
          p_orth_space_basis->push_back(constr_vec(
            m_ambient_dim,
            eig.eigenvectors().col(j).data(),
            eig.eigenvectors().col(j).data() + m_ambient_dim));
        }
      }
    }

    m_are_tangent_spaces_computed[i] = true;

    //*************************************************************************

    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, FT(1)/sqrt(sqlen(n)));
    //std::cerr << "IP = " << inner_pdct(n, ts[0]) << " & " << inner_pdct(n, ts[1]) << std::endl;

    return tsb;
    
#endif

    /*
    // Alternative code (to be used later)
    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, FT(1)/sqrt(sqlen(n)));
    //Vector t1(12., 15., 65.);
    //Vector t2(32., 5., 85.);
    //Tangent_space_basis ts;
    //ts.reserve(m_intrinsic_dimension);
    //ts.push_back(diff_vec(t1, scaled_vec(n, inner_pdct(t1, n))));
    //ts.push_back(diff_vec(t2, scaled_vec(n, inner_pdct(t2, n))));
    //ts = compute_gram_schmidt_basis(ts, m_k);
    //return ts;
    */
  }

  Point compute_perturbed_point(std::size_t pt_idx) const
  {
#ifdef CGAL_TC_PERTURB_POSITION
    return m_k.translated_point_d_object()(
      m_points[pt_idx], m_translations[pt_idx]);
#else
    return m_points[pt_idx];
#endif
  }

  void compute_perturbed_weighted_point(std::size_t pt_idx, Point &p, FT &w) const
  {
#ifdef CGAL_TC_PERTURB_POSITION
    p = m_k.translated_point_d_object()(
      m_points[pt_idx], m_translations[pt_idx]);
#else
    p = m_points[pt_idx];
#endif
    w = m_weights[pt_idx];
  }

  Weighted_point compute_perturbed_weighted_point(std::size_t pt_idx) const
  {
    typename Kernel::Construct_weighted_point_d k_constr_wp =
      m_k.construct_weighted_point_d_object();

    Weighted_point wp = k_constr_wp(
#ifdef CGAL_TC_PERTURB_POSITION
      m_k.translated_point_d_object()(m_points[pt_idx], m_translations[pt_idx]),
#else
      m_points[pt_idx],
#endif
      m_weights[pt_idx]);

    return wp;
  }

  Point unproject_point(const Tr_point &p,
                        const Tangent_space_basis &tsb,
                        const Tr_traits &tr_traits) const
  {
    typename Kernel::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename Kernel::Scaled_vector_d k_scaled_vec =
      m_k.scaled_vector_d_object();
    typename Tr_traits::Compute_coordinate_d coord =
      tr_traits.compute_coordinate_d_object();

    Point global_point = tsb.origin();
    for (int i = 0 ; i < m_intrinsic_dimension ; ++i)
      global_point = k_transl(global_point,
                              k_scaled_vec(tsb[i], coord(p, i)));

    return global_point;
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_bare_point project_point(const Point &p,
                              const Tangent_space_basis &tsb) const
  {
    typename Kernel::Scalar_product_d inner_pdct =
      m_k.scalar_product_d_object();
    typename Kernel::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    coords.reserve(m_intrinsic_dimension);
    for (std::size_t i = 0 ; i < m_intrinsic_dimension ; ++i)
    {
      // Local coords are given by the inner product with the vectors of tsb
      Vector v = diff_points(p, tsb.origin());
      FT coord = inner_pdct(v, tsb[i]);
      coords.push_back(coord);
    }

    return Tr_bare_point(m_intrinsic_dimension, coords.begin(), coords.end());
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_point project_point_and_compute_weight(const Weighted_point &wp,
                                            const Tangent_space_basis &tsb,
                                            const Tr_traits &tr_traits) const
  {
    typename Kernel::Point_drop_weight_d k_drop_w =
      m_k.point_drop_weight_d_object();
    typename Kernel::Point_weight_d k_point_weight =
      m_k.point_weight_d_object();
    return project_point_and_compute_weight(
      k_drop_w(wp), k_point_weight(wp), tsb, tr_traits);
  }

  Tr_point project_point_and_compute_weight(const Point &p, const FT w,
                                            const Tangent_space_basis &tsb,
                                            const Tr_traits &tr_traits) const
  {
    const int point_dim = m_k.point_dimension_d_object()(p);
    typename Kernel::Scalar_product_d inner_pdct =
      m_k.scalar_product_d_object();
    typename Kernel::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();
    typename Kernel::Construct_cartesian_const_iterator_d ccci =
      m_k.construct_cartesian_const_iterator_d_object();

    Vector v = diff_points(p, tsb.origin());

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(ccci(tsb.origin()), ccci(tsb.origin(), 0));
    coords.reserve(m_intrinsic_dimension);
    for (std::size_t i = 0 ; i < m_intrinsic_dimension ; ++i)
    {
      // Local coords are given by the inner product with the vectors of tsb
      FT coord = inner_pdct(v, tsb[i]);
      coords.push_back(coord);

      // p_proj += coord * v;
      for (int j = 0 ; j < point_dim ; ++j)
        p_proj[j] += coord * tsb[i][j];
    }

    Point projected_pt(point_dim, p_proj.begin(), p_proj.end());

    return tr_traits.construct_weighted_point_d_object()
    (
      tr_traits.construct_point_d_object()(
        m_intrinsic_dimension, coords.begin(), coords.end()),
      w - m_k.squared_distance_d_object()(p, projected_pt)
    );
  }

  // A simplex here is a local tri's full cell handle
  bool is_simplex_consistent(Tr_full_cell_handle fch, int cur_dim) const
  {
    std::set<std::size_t> c;
    for (int i = 0 ; i < cur_dim + 1 ; ++i)
    {
      std::size_t data = fch->vertex(i)->data();
      c.insert(data);
    }
    return is_simplex_consistent(c);
  }

  // A simplex here is a list of point indices
  bool is_simplex_consistent(std::set<std::size_t> const& simplex) const
  {
    int cur_dim_plus_1 = static_cast<int>(simplex.size());

    // Check if the simplex is in the stars of all its vertices
    std::set<std::size_t>::const_iterator it_point_idx = simplex.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for ( ; it_point_idx != simplex.end() ; ++it_point_idx)
    {
      std::size_t point_idx = *it_point_idx;
      // Don't check infinite simplices
      if (point_idx == std::numeric_limits<std::size_t>::max())
        continue;

      Star const& star = m_stars[point_idx];

      // What we're looking for is "simplex" \ point_idx
      Incident_simplex ic_to_find = simplex;
      ic_to_find.erase(point_idx);

      // For each cell
      if (std::find(star.begin(), star.end(), ic_to_find) == star.end())
        return false;
    }

    return true;
  }

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for try_to_solve_inconsistencies_in_a_local_triangulation function
  class Try_to_solve_inconsistencies_in_a_local_triangulation
  {
    Tangential_complex & m_tc;
    tbb::combinable<std::size_t> &m_num_inconsistencies;

  public:
    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
      Tangential_complex &tc, tbb::combinable<std::size_t> &num_inconsistencies)
    : m_tc(tc), m_num_inconsistencies(num_inconsistencies)
    {}

    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
      const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc), m_num_inconsistencies(ctt.m_num_inc)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
      {
        m_num_inconsistencies.local() +=
          m_tc.try_to_solve_inconsistencies_in_a_local_triangulation(i);
      }
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  void perturb(std::size_t point_idx)
  {
    // Perturb the weight?
#ifdef CGAL_TC_PERTURB_WEIGHT
    m_weights[point_idx] = m_random_generator.get_double(0., m_sq_half_sparsity);
    if(m_weights_memory.size() > 0) // external weights were initially set
      m_weights[point_idx] = m_weights[point_idx] + m_weights_memory[point_idx];
#endif

#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    m_perturb_tangent_space[point_idx] = true;
#endif

    // Perturb the position?
#ifdef CGAL_TC_PERTURB_POSITION
# ifdef CGAL_TC_PERTURB_POSITION_GLOBAL
    typename Kernel::Point_to_vector_d k_pt_to_vec =
      m_k.point_to_vector_d_object();
    CGAL::Random_points_on_sphere_d<Point>
      tr_point_on_sphere_generator(
        m_ambient_dim, m_random_generator.get_double(0., m_half_sparsity));
    // Parallel
#  if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
    Vector transl = k_pt_to_vec(*tr_point_on_sphere_generator);
    m_p_perturb_mutexes[point_idx].lock();
    m_translations[point_idx] = transl;
    m_p_perturb_mutexes[point_idx].unlock();
    // Sequential
#  else
    m_translations[point_idx] = k_pt_to_vec(*tr_point_on_sphere_generator);
#  endif

# else // CGAL_TC_PERTURB_POSITION_TANGENTIAL
    const Tr_traits &local_tr_traits =
      m_triangulations[point_idx].tr().geom_traits();
    typename Tr_traits::Compute_coordinate_d coord =
      local_tr_traits.compute_coordinate_d_object();
    typename Kernel::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename Kernel::Construct_vector_d k_constr_vec =
      m_k.construct_vector_d_object();
    typename Kernel::Scaled_vector_d k_scaled_vec =
      m_k.scaled_vector_d_object();

    CGAL::Random_points_on_sphere_d<Tr_bare_point> 
      tr_point_on_sphere_generator(
        m_intrinsic_dimension, 
        m_random_generator.get_double(0., m_half_sparsity));

    Tr_point local_random_transl =
      local_tr_traits.construct_weighted_point_d_object()(
        *tr_point_on_sphere_generator, 0);
    Translation_for_perturb global_transl = k_constr_vec(m_ambient_dim);
    const Tangent_space_basis &tsb = m_tangent_spaces[point_idx];
    for (int i = 0 ; i < m_intrinsic_dimension ; ++i)
    {
      global_transl = k_transl(
        global_transl,
        k_scaled_vec(tsb[i], coord(local_random_transl, i))
      );
    }
    // Parallel
#  if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
    m_p_perturb_mutexes[point_idx].lock();
    m_translations[point_idx] = global_transl;
    m_p_perturb_mutexes[point_idx].unlock();
    // Sequential
#  else
    m_translations[point_idx] = global_transl;
#  endif

# endif // CGAL_TC_PERTURB_POSITION_TANGENTIAL
#endif // CGAL_TC_PERTURB_POSITION
  }

  bool try_to_solve_inconsistencies_in_a_local_triangulation(
                                                          std::size_t tr_index)
  {
    bool is_inconsistent = false;

#ifdef CGAL_LINKED_WITH_TBB
    //Tr_mutex::scoped_lock lock(m_tr_mutexes[tr_index]);
#endif

    Star const& star = m_stars[tr_index];
    Triangulation const& tr    = m_triangulations[tr_index].tr();
    Tr_vertex_handle center_vh = m_triangulations[tr_index].center_vertex();
    const Tr_traits &local_tr_traits = tr.geom_traits();
    int cur_dim = tr.current_dimension();

    // For each incident simplex
    Star::const_iterator it_inc_simplex = star.begin();
    Star::const_iterator it_inc_simplex_end = star.end();
    for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
    {
      const Incident_simplex &incident_simplex = *it_inc_simplex;

      // Don't check infinite cells
      if (*incident_simplex.rbegin() == std::numeric_limits<std::size_t>::max())
        continue;

      std::set<std::size_t> c = incident_simplex;
      c.insert(tr_index); // Add the missing index

//*****************************************************************************
// STRATEGY 1: perturb all the points of the first inconsistent simplex
//*****************************************************************************
#ifdef CGAL_TC_PERTURB_THE_SIMPLEX_ONLY
      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        for (std::set<std::size_t>::const_iterator it = c.begin();
             it != c.end() ; ++it)
        {
          perturb(*it);
        }

# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time
        break;
      }

//*****************************************************************************
// STRATEGY 2: perturb the center point only
//*****************************************************************************
#elif defined(CGAL_TC_PERTURB_THE_CENTER_VERTEX_ONLY)
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        std::size_t idx = tr_index;
        /*int k;
        do
        {
          k = rand() % tr.current_dimension();
        } while ((*it_c)->vertex(k) == center_vh);
        std::size_t idx = (*it_c)->vertex(k)->data();*/

        perturb(idx);

# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time
        break;
      }

//*****************************************************************************
// STRATEGY 3: perturb all the points of the 1-star
//*****************************************************************************
#elif defined(CGAL_TC_PERTURB_THE_1_STAR)

      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        std::set<std::size_t> the_1_star;

        Star::const_iterator it_inc_simplex = star.begin();
        Star::const_iterator it_inc_simplex_end = star.end();
        for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
        {
          the_1_star.insert(it_inc_simplex->begin(), it_inc_simplex ->end());
        }

        for (std::set<std::size_t>::iterator it = the_1_star.begin() ;
             it != the_1_star.end() ; ++it)
        {
          perturb(*it);
        }

# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time
        break;
      }

//*****************************************************************************
// STRATEGY 4: perturb the k + 1 + CGAL_TC_NUMBER_OF_ADDITIONNAL_PERTURBED_POINTS
// closest points (to the power center of first the inconsistent cell)
//*****************************************************************************
#elif defined(CGAL_TC_PERTURB_N_CLOSEST_POINTS)

      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        // Get the k + 1 + CGAL_TC_NUMBER_OF_ADDITIONNAL_PERTURBED_POINTS
        // closest points

        std::vector<Tr_point> simplex_pts;
        simplex_pts.reserve(c.size());

        Incident_simplex::const_iterator it_point_idx = c.begin();
        Incident_simplex::const_iterator it_point_idx_end = c.end();
        // For each point p of the simplex, we reproject it onto the tangent
        // space. Could be optimized since it's already been computed before.
        for ( ; it_point_idx != it_point_idx_end ; ++it_point_idx)
        {
          simplex_pts.push_back(project_point_and_compute_weight(
            m_points[*it_point_idx], m_weights[*it_point_idx],
              m_tangent_spaces[tr_index], local_tr_traits));
        }

        typename Tr_traits::Power_center_d power_center =
          local_tr_traits.power_center_d_object();
        typename Tr_traits::Compute_coordinate_d coord =
          local_tr_traits.compute_coordinate_d_object();

        Point global_center = unproject_point(
          power_center(simplex_pts.begin(), simplex_pts.end()),
          m_tangent_spaces[tr_index],
          local_tr_traits);

        KNS_range kns_range = m_points_ds.query_ANN(
          global_center,
          CGAL_TC_NUMBER_OF_PERTURBED_POINTS(m_intrinsic_dimension));
        std::vector<std::size_t> neighbors;
        for (KNS_iterator nn_it = kns_range.begin() ;
             nn_it != kns_range.end() ;
             ++nn_it)
        {
          neighbors.push_back(nn_it->first);
        }

        for (std::vector<std::size_t>::iterator it = neighbors.begin();
             it != neighbors.end() ;
             ++it)
        {
          perturb(*it);
        }

# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time
        break;
      }
//*****************************************************************************
// STRATEGY 5: perturb one random point of the simplex
//*****************************************************************************
#else
      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;
        int rnd = m_random_generator.get_int(0, static_cast<int>(c.size()));
        if (rnd == 0)
          perturb(tr_index);
        else
        {
          std::set<std::size_t>::const_iterator it_idx = c.begin();
          std::advance(it_idx, rnd - 1);
          perturb(*it_idx);
        }

# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time
        break;
      }

#endif // CGAL_TC_PERTURB_THE_SIMPLEX_ONLY
    }

    return is_inconsistent;
  }

  std::ostream &export_vertices_to_off(
    std::ostream & os, std::size_t &num_vertices,
    bool use_perturbed_points = false) const
  {
    if (m_points.empty())
    {
      num_vertices = 0;
      return os;
    }

    // If m_intrinsic_dimension = 1, we output each point two times
    // to be able to export each segment as a flat triangle with 3 different
    // indices (otherwise, Meshlab detects degenerated simplices)
    const int N = (m_intrinsic_dimension == 1 ? 2 : 1);

    // Kernel functors
    typename Kernel::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();

    int num_coords = min(m_ambient_dim, 3);
#ifdef CGAL_TC_EXPORT_NORMALS
    OS_container::const_iterator it_os = m_orth_spaces.begin();
#endif
    typename Points::const_iterator it_p = m_points.begin();
    typename Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for (std::size_t i = 0 ; it_p != it_p_end ; ++it_p, ++i)
    {
      Point p = (use_perturbed_points ? compute_perturbed_point(i) : *it_p);
      for (int ii = 0 ; ii < N ; ++ii)
      {
        int i = 0;
        for ( ; i < num_coords ; ++i)
          os << CGAL::to_double(coord(p, i)) << " ";
        if (i == 2)
          os << "0";

#ifdef CGAL_TC_EXPORT_NORMALS
        for (i = 0 ; i < num_coords ; ++i)
          os << " " << CGAL::to_double(coord(*it_os->begin(), i));
#endif
        os << std::endl;
      }
#ifdef CGAL_TC_EXPORT_NORMALS
      ++it_os;
#endif
    }

    num_vertices = N*m_points.size();
    return os;
  }

  void insert_higher_dim_simplex_into_star(
    std::size_t index,
    const std::set<std::size_t> &simplex)
  {
    Incident_simplex incident_simplex = simplex;
    incident_simplex.erase(index); // Remove the center index

    Star &star = m_stars[index];

    std::set<std::size_t>::const_iterator it_point_idx = simplex.begin();
    std::set<std::size_t>::const_iterator it_point_idx_end = simplex.end();
    for ( ; it_point_idx != it_point_idx_end ; ++it_point_idx)
    {
      // Skip center index
      if (*it_point_idx == index)
        continue;

      // Temporarily remove this index
      incident_simplex.erase(*it_point_idx);
      // Erase incident_simplex from star
      star.erase(std::remove(star.begin(), star.end(), incident_simplex),
                 star.end());
      incident_simplex.insert(*it_point_idx);
    }

    star.push_back(incident_simplex);
  }

  // Solves one inconsistency
  // "inconsistent_simplex" must contain p_idx and q_idx
  // "inconsistent_simplex" must be in star(p) but not in star(q)
  void solve_inconsistency_by_adding_higher_dimensional_simplices(
    std::size_t p_idx, std::size_t q_idx,
    const std::set<std::size_t> &inconsistent_simplex)
  {
    CGAL_assertion_code(
      std::set<std::size_t> inc_s_minus_p = inconsistent_simplex;
      inc_s_minus_p.erase(p_idx);
      std::set<std::size_t> inc_s_minus_q = inconsistent_simplex;
      inc_s_minus_q.erase(q_idx);
    );
    CGAL_assertion(std::find(m_stars[p_idx].begin(), m_stars[p_idx].end(),
                             inc_s_minus_p) != m_stars[p_idx].end());
    CGAL_assertion(std::find(m_stars[q_idx].begin(), m_stars[q_idx].end(),
                             inc_s_minus_q) == m_stars[q_idx].end());

    typename Kernel::Point_drop_weight_d k_drop_w =
      m_k.point_drop_weight_d_object();
    typename Kernel::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename Kernel::Squared_distance_d k_sqdist =
      m_k.squared_distance_d_object();
    typename Kernel::Difference_of_points_d k_diff_pts =
      m_k.difference_of_points_d_object();
    typename Kernel::Scalar_product_d k_inner_pdct =
      m_k.scalar_product_d_object();
    typename Kernel::Construct_weighted_point_d k_constr_wp =
      m_k.construct_weighted_point_d_object();
    typename Kernel::Power_distance_d k_power_dist =
      m_k.power_distance_d_object();

    const Tr_traits &q_tr_traits = m_triangulations[q_idx].tr().geom_traits();
    typename Tr_traits::Power_center_d tr_power_center =
      q_tr_traits.power_center_d_object();
    typename Tr_traits::Point_weight_d tr_point_weight =
      q_tr_traits.point_weight_d_object();

    //-------------------------------------------------------------------------
    //1. Compute power_center(p'q'r1'r2'..ri') in Tp => Cp
    //2. Compute power_center(inconsistent_simplex projected in Tq)
    // => gives Cq and radius Rq
    // Rq is also the radius of the ambient sphere S whose center is Cq and
    // which goes through all the ambient points of "inconsistent_simplex"
    //------------------------------------------------------------------------
    std::vector<Tr_point> simplex_pts_in_Tp, simplex_pts_in_Tq;
    simplex_pts_in_Tp.reserve(inconsistent_simplex.size());
    simplex_pts_in_Tq.reserve(inconsistent_simplex.size());

    // No need to lock the mutex here since this will not be called while
    // other threads are perturbing the positions
    const Point pt_p = compute_perturbed_point(p_idx);

    std::set<std::size_t>::const_iterator it_point_idx =
                                                  inconsistent_simplex.begin();
    std::set<std::size_t>::const_iterator it_point_idx_end =
                                                  inconsistent_simplex.end();
    // For each point of the simplex, we reproject it onto the tangent
    // space. Could be optimized since it's already been computed before.
    for ( ; it_point_idx != it_point_idx_end ; ++it_point_idx)
    {
      const Weighted_point wp = compute_perturbed_weighted_point(*it_point_idx);
      // No need to lock the Mutex_for_perturb here since this will not be
      // called while other threads are perturbing the positions
      simplex_pts_in_Tp.push_back(project_point_and_compute_weight(
        wp, m_tangent_spaces[p_idx], q_tr_traits));
      simplex_pts_in_Tq.push_back(project_point_and_compute_weight(
        wp, m_tangent_spaces[q_idx], q_tr_traits));
    }

    Tr_point Cp = tr_power_center(
      simplex_pts_in_Tp.begin(), simplex_pts_in_Tp.end());
    Tr_point Cq = tr_power_center(
      simplex_pts_in_Tq.begin(), simplex_pts_in_Tq.end());

    FT circumsphere_sqradius_p = tr_point_weight(Cp);
    FT circumsphere_sqradius_q = tr_point_weight(Cq);
#ifdef CGAL_TC_PERTURB_WEIGHT
    FT squared_circumsphere_radius_q_plus_margin =
      circumsphere_sqradius_q + 4*m_sq_half_sparsity;
#else
    FT squared_circumsphere_radius_q_plus_margin = CGAL::square(
      CGAL::sqrt(circumsphere_sqradius_q) + 2*m_half_sparsity);
#endif

    Weighted_point global_Cp = k_constr_wp(
      unproject_point(Cp, m_tangent_spaces[p_idx], q_tr_traits),
      circumsphere_sqradius_p);

    Weighted_point global_Cq = k_constr_wp(
      unproject_point(Cq, m_tangent_spaces[q_idx], q_tr_traits),
      circumsphere_sqradius_q);

    // CJTODO TEMP ====================
    /*{
    INS_range ins_range = m_points_ds.query_incremental_ANN(k_drop_w(global_Cp));
    for (INS_iterator nn_it = ins_range.begin() ;
         nn_it != ins_range.end() ;
         ++nn_it)
    {
      FT neighbor_sqdist = nn_it->second;

      //std::cerr << nn_it->first << " : " << neighbor_sqdist << " / "; // CJTODO TEMP

      // When we're sure we got all the potential points, break
      //if (neighbor_sqdist > circumsphere_sqradius_p + m_sq_half_sparsity)
      //  break;

      std::size_t neighbor_point_idx = nn_it->first;
      FT point_to_Cp_power_sqdist = k_power_dist(
        global_Cp, compute_perturbed_weighted_point(neighbor_point_idx));
      //std::cerr << point_to_Cp_power_sqdist << std::endl; // CJTODO TEMP
      // If the point is ACTUALLY "inside" S
      if (point_to_Cp_power_sqdist <= FT(0)
        && inconsistent_simplex.find(neighbor_point_idx) ==
           inconsistent_simplex.end())
      {
        std::cerr << "Warning: " << neighbor_point_idx << " is inside Cp with power dist " << point_to_Cp_power_sqdist << "\n";
      }
    }
    }*/
    // /CJTODO ====================

    //-------------------------------------------------------------------------
    //3. Find points t1, t2... (in ambient space) which are inside S
    //-------------------------------------------------------------------------
    std::vector<std::size_t> inside_pt_indices;
    INS_range ins_range = m_points_ds.query_incremental_ANN(k_drop_w(global_Cq));
    for (INS_iterator nn_it = ins_range.begin() ;
         nn_it != ins_range.end() ;
         ++nn_it)
    {
      FT neighbor_sqdist = nn_it->second;

      // When we're sure we got all the potential points, break
      if (neighbor_sqdist > squared_circumsphere_radius_q_plus_margin)
        break;

      std::size_t neighbor_point_idx = nn_it->first;
      FT point_to_Cq_power_sqdist = k_power_dist(
        global_Cq, compute_perturbed_weighted_point(neighbor_point_idx));
      // If the point is ACTUALLY "inside" S
      if (point_to_Cq_power_sqdist <= FT(0)
        && inconsistent_simplex.find(neighbor_point_idx) ==
           inconsistent_simplex.end())
      {
        inside_pt_indices.push_back(neighbor_point_idx);
      }
      // CJTODO: use this instead of point_to_Cq_power_sqdist?
      /*{
        typename Tr_traits::Power_test_d side = q_tr_traits.power_test_d_object();
        typename Tr_traits::Orientation_d orient = q_tr_traits.orientation_d_object();
        Orientation o = orient(simplex_pts_in_Tq.begin(), simplex_pts_in_Tq.end());
        auto p = project_point_and_compute_weight(
          compute_perturbed_weighted_point(neighbor_point_idx),
          m_tangent_spaces[q_idx], q_tr_traits);
        auto sid = (o == NEGATIVE ?
          side(simplex_pts_in_Tq.rbegin(), simplex_pts_in_Tq.rend(), p)
          : side(simplex_pts_in_Tq.begin(), simplex_pts_in_Tq.end(), p));
        switch(sid)
        {
        case ON_NEGATIVE_SIDE:
          std::cerr << "ON_NEGATIVE_SIDE" << std::endl; // CJTODO TEMP
          break;
        case ON_POSITIVE_SIDE:
          std::cerr << "ON_POSITIVE_SIDE" << std::endl; // CJTODO TEMP
          break;
        case ON_ORIENTED_BOUNDARY:
          std::cerr << "ON_ORIENTED_BOUNDARY" << std::endl; // CJTODO TEMP
          break;
        }
      }*/
    }
    CGAL_assertion_msg(!inside_pt_indices.empty(),
      "There should be at least one vertex inside the sphere");

    // CJTODO TEMP DEBUG
    /*if (inside_pt_indices.empty())
    {
      //compute_tangent_triangulation(q_idx, true);
      std::cerr << "Error: inside_pt_indices.empty()\n";
      std::cerr << "Stars:\n";
      for (auto s : m_stars[q_idx])
      {
        std::cerr << q_idx << " ";
        std::copy(s.begin(), s.end(),
          std::ostream_iterator<std::size_t>(std::cerr, " "));
        std::cerr << std::endl;
      }
      std::cerr << std::endl;
    }*/
    // CJTODO TEMP DEBUG
    if (inside_pt_indices.size() > 1)
    {
      std::cerr << "Warning: " << inside_pt_indices.size() << " insiders in "
        << inconsistent_simplex.size() - 1 << " simplex\n";
    }

    //-------------------------------------------------------------------------
    //4. If there's more than one ti... or not
    //-------------------------------------------------------------------------
    std::size_t inside_point_idx;
    if (inside_pt_indices.size() > 1)
    {
      //-----------------------------------------------------------------------
      //5. For each ti, compute the sphere that goes through
      //   p, q, r1, r2..., ri and ti whose center is on (cp, cq)
      //   We're looking for a point on (Cp, Cq) at equal distance from p and
      //   ti.
      //   The center of the sphere is then: Cp + a(Cq - Cp)
      //   where a = (sqdist(Cp,ti) - sqdist(Cp,p)) / (2*(Cq-Cp).(ti-p))
      //6. Keep point ti such as dist(cp, ci) is the smallest
      //-----------------------------------------------------------------------

      FT min_a = std::numeric_limits<FT>::max();
      for (std::size_t i = 0 ; i < inside_pt_indices.size() ; ++i)
      {
        std::size_t idx = inside_pt_indices[i];
        const Point ti = compute_perturbed_point(idx);
        const Point &cp = k_drop_w(global_Cp);
        const Point &cq = k_drop_w(global_Cq);

#ifdef CGAL_TC_PERTURB_WEIGHT
        const Weighted_point ti_w = compute_perturbed_weighted_point(idx);
        const Weighted_point p_w = compute_perturbed_weighted_point(p_idx);
        const Weighted_point cp_w0 = k_constr_wp(k_drop_w(global_Cp), FT(0));
        const Weighted_point wp_w0 = k_constr_wp(k_drop_w(global_Cq), FT(0));
        FT a =
          (k_power_dist(cp_w0, ti_w) - k_power_dist(cp_w0, p_w)) /
          (FT(2)*k_inner_pdct(k_diff_pts(cq, cp), k_diff_pts(ti, pt_p)));
#else
        FT a =
          (k_sqdist(cp, ti) - k_sqdist(cp, pt_p)) /
          (FT(2)*k_inner_pdct(k_diff_pts(cq, cp), k_diff_pts(ti, pt_p)));
#endif

        if (a < min_a)
        {
          min_a = a;
          inside_point_idx = idx;
        }
      }

      // CJTODO TEMP ====================
      /*{
      typename Kernel::Scaled_vector_d scaled_vec = m_k.scaled_vector_d_object();
      typename Kernel::Point_weight_d k_weight = m_k.point_weight_d_object();
      Weighted_point C = k_constr_wp(
        k_transl(k_drop_w(global_Cp), scaled_vec(k_diff_pts(k_drop_w(global_Cq), k_drop_w(global_Cp)), min_a)),
        k_sqdist(k_transl(k_drop_w(global_Cp), scaled_vec(k_diff_pts(k_drop_w(global_Cq), k_drop_w(global_Cp)), min_a)), pt_p));
      INS_range ins_range = m_points_ds.query_incremental_ANN(k_drop_w(C));
      for (INS_iterator nn_it = ins_range.begin() ;
           nn_it != ins_range.end() ;
           ++nn_it)
      {
        FT neighbor_sqdist = nn_it->second;

        //std::cerr << nn_it->first << " : " << neighbor_sqdist << " / "; // CJTODO TEMP

        // When we're sure we got all the potential points, break
        if (neighbor_sqdist > k_weight(C) + m_sq_half_sparsity)
          break;

        std::size_t neighbor_point_idx = nn_it->first;
        FT point_to_C_power_sqdist =
          k_power_dist(C, compute_perturbed_weighted_point(neighbor_point_idx));
        //std::cerr << point_to_Cp_power_sqdist << std::endl; // CJTODO TEMP
        // If the point is ACTUALLY "inside" S
        if (point_to_C_power_sqdist <= FT(-0.000001)
          && inconsistent_simplex.find(neighbor_point_idx) ==
             inconsistent_simplex.end())
        {
          std::cerr << "Warning: " << neighbor_point_idx << " is inside C with power dist " << point_to_C_power_sqdist << "\n";
        }
      }
      }*/
      // /CJTODO ====================
    }
    else
    {
      inside_point_idx = *inside_pt_indices.begin();
    }

    //-------------------------------------------------------------------------
    //7. Create a k+1-simplex (inconsistent_simplex, ti)
    //-------------------------------------------------------------------------
    std::set<std::size_t> new_simplex = inconsistent_simplex;
    new_simplex.insert(inside_point_idx);

    it_point_idx = new_simplex.begin();
    it_point_idx_end = new_simplex.end();
    for ( ; it_point_idx != it_point_idx_end ; ++it_point_idx)
    {
      insert_higher_dim_simplex_into_star(*it_point_idx, new_simplex);
    }
    // CJTODO: call
    // check_and_solve_inconsistencies_by_adding_higher_dim_simplices
    // recursively? Not sure, since the star will be parsed again from
    // the beginning
  }

  // Test and solve inconsistencies of a simplex.
  // Returns true if some inconsistencies were found.
  // Precondition: incident_simplex is in the star of m_points[tr_index]
  bool check_and_solve_inconsistencies_by_adding_higher_dim_simplices(
    std::size_t tr_index, const std::set<std::size_t> &incident_simplex)
  {
    bool inconsistencies_found = false;

    // Don't check infinite simplices
    if (*incident_simplex.rbegin()
        == std::numeric_limits<std::size_t>::max())
      return false;

    std::set<std::size_t> simplex = incident_simplex;
    simplex.insert(tr_index);

    // Check if the simplex is in the stars of all its vertices
    std::set<std::size_t>::const_iterator it_point_idx = incident_simplex.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for ( ; it_point_idx != incident_simplex.end() ; ++it_point_idx)
    {
      std::size_t point_idx = *it_point_idx;

      Star const& star = m_stars[point_idx];

      // What we're looking for is "simplex" \ point_idx
      Incident_simplex ic_to_find = simplex;
      ic_to_find.erase(point_idx);

      if (std::find(star.begin(), star.end(), ic_to_find) == star.end())
      {
        solve_inconsistency_by_adding_higher_dimensional_simplices(
          tr_index, *it_point_idx, simplex);
        inconsistencies_found = true;
        break;
      }
    }

    return inconsistencies_found;
  }

  std::ostream &export_simplices_to_off(
    std::ostream & os, std::size_t &num_simplices,
    bool color_inconsistencies = false,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_red = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_green = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_blue = NULL)
    const
  {
    // If m_intrinsic_dimension = 1, each point is output two times
    // (see export_vertices_to_off)
    num_simplices = 0;
    std::size_t num_inconsistent_simplices = 0;
    std::size_t num_inconsistent_stars = 0;
    typename Tr_container::const_iterator it_tr = m_triangulations.begin();
    typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for (std::size_t idx = 0 ; it_tr != it_tr_end ; ++it_tr, ++idx)
    {
      bool is_star_inconsistent = false;

      Triangulation const& tr    = it_tr->tr();
      Tr_vertex_handle center_vh = it_tr->center_vertex();

      if (&tr == NULL || tr.current_dimension() < m_intrinsic_dimension)
        continue;

      // Color for this star
      std::stringstream color;
      //color << rand()%256 << " " << 100+rand()%156 << " " << 100+rand()%156;
      color << 128 << " " << 128 << " " << 128;

      // Gather the triangles here, with an int telling its color
      typedef std::vector<std::pair<std::set<std::size_t>, int> >
                                                           Star_using_triangles;
      Star_using_triangles star_using_triangles;

      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        std::set<std::size_t> c = *it_inc_simplex;
        c.insert(idx);
        std::size_t num_vertices = c.size();

        int color_simplex = -1;// -1=no color, 0=yellow, 1=red, 2=green, 3=blue
        if (color_inconsistencies)
        {
          is_star_inconsistent = !is_simplex_consistent(c);
          color_simplex = (is_star_inconsistent ? 0 : -1);
        }

        if (color_simplex == -1)
        {
          if (p_simpl_to_color_in_red && 
              std::find(
                p_simpl_to_color_in_red->begin(),
                p_simpl_to_color_in_red->end(),
                c) != p_simpl_to_color_in_red->end())
          {
            color_simplex = 1;
          }
          else if (p_simpl_to_color_in_green && 
              std::find(
                p_simpl_to_color_in_green->begin(),
                p_simpl_to_color_in_green->end(),
                c) != p_simpl_to_color_in_green->end())
          {
            color_simplex = 2;
          }
          else if (p_simpl_to_color_in_blue && 
              std::find(
                p_simpl_to_color_in_blue->begin(),
                p_simpl_to_color_in_blue->end(),
                c) != p_simpl_to_color_in_blue->end())
          {
            color_simplex = 3;
          }
        }

        // If m_intrinsic_dimension = 1, each point is output two times,
        // so we need to multiply each index by 2
        // And if only 2 vertices, add a third one (each vertex is duplicated in
        // the file when m_intrinsic dim = 2)
        if (m_intrinsic_dimension == 1)
        {
          std::set<std::size_t> tmp_c;
          std::set<std::size_t>::iterator it = c.begin();
          for ( ; it != c.end() ; ++it)
            tmp_c.insert(*it * 2);
          if (num_vertices == 2)
            tmp_c.insert(*tmp_c.rbegin() + 1);

          c = tmp_c;
        }

        if (num_vertices <= 3)
        {
          star_using_triangles.push_back(std::make_pair(c, color_simplex));
        }
        else
        {
          // num_vertices >= 4: decompose the simplex in triangles
          std::vector<bool> booleans(num_vertices, false);
          std::fill(booleans.begin() + num_vertices - 3, booleans.end(), true);
          do
          {
            std::set<std::size_t> triangle;
            std::set<std::size_t>::iterator it = c.begin();
            for (int i = 0; it != c.end() ; ++i, ++it)
            {
              if (booleans[i])
                triangle.insert(*it);
            }
            star_using_triangles.push_back(
              std::make_pair(triangle, color_simplex));
          } while (std::next_permutation(booleans.begin(), booleans.end()));
        }
      }

      // For each cell
      Star_using_triangles::const_iterator it_simplex =
                                                  star_using_triangles.begin();
      Star_using_triangles::const_iterator it_simplex_end =
                                                  star_using_triangles.end();
      for ( ; it_simplex != it_simplex_end ; ++it_simplex)
      {
        // Don't export infinite cells
        if (*it_simplex->first.rbegin()
            == std::numeric_limits<std::size_t>::max())
          continue;

        const std::set<std::size_t> &c = it_simplex->first;
        int color_simplex = it_simplex->second;

        std::stringstream sstr_c;

        std::set<std::size_t>::const_iterator it_point_idx = c.begin();
        for ( ; it_point_idx != c.end() ; ++it_point_idx)
        {
          sstr_c << *it_point_idx << " ";
        }

        // In order to have only one time each simplex, we only keep it
        // if the lowest index is the index of the center vertex
        if (*c.begin() != idx && color_simplex == -1)
          continue;

        os << 3 << " " << sstr_c.str();
        if (color_inconsistencies || p_simpl_to_color_in_red
            || p_simpl_to_color_in_green || p_simpl_to_color_in_blue)
        {
          switch (color_simplex)
          {
            case 0: os << " 255 255 0"; ++num_inconsistent_simplices; break;
            case 1: os << " 255 0 0"; break;
            case 2: os << " 0 255 0"; break;
            case 3: os << " 0 0 255"; break;
            default: os << " " << color.str(); break;
          }            
        }
        ++num_simplices;
        os << std::endl;
      }
      if (is_star_inconsistent)
        ++num_inconsistent_stars;
    }

#ifdef CGAL_TC_VERBOSE
    std::cerr << std::endl
      << "=========================================================="
      << std::endl
      << "Export to OFF:\n"
      << "  * Number of vertices: " << m_points.size() << std::endl
      << "  * Total number of simplices: " << num_simplices << std::endl
      << "  * Number of inconsistent stars: "
      << num_inconsistent_stars << " ("
      << (m_points.size() > 0 ?
          100. * num_inconsistent_stars / m_points.size() : 0.) << "%)"
      << std::endl
      << "  * Number of inconsistent simplices: "
      << num_inconsistent_simplices << " ("
      << (num_simplices > 0 ?
          100. * num_inconsistent_simplices / num_simplices : 0.) << "%)"
      << std::endl
      << "=========================================================="
      << std::endl;
#endif

    return os;
  }

public:
  std::ostream &export_simplices_to_off(
    const Simplicial_complex &complex,
    std::ostream & os, std::size_t &num_simplices,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_red = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_green = NULL,
    std::set<std::set<std::size_t> > const *p_simpl_to_color_in_blue = NULL)
    const
  {
    typedef Simplicial_complex::Simplex                     Simplex;
    typedef Simplicial_complex::Simplex_range               Simplex_range;

    // If m_intrinsic_dimension = 1, each point is output two times
    // (see export_vertices_to_off)
    num_simplices = 0;
    typename Simplex_range::const_iterator it_s =
      complex.simplex_range().begin();
    typename Simplex_range::const_iterator it_s_end =
      complex.simplex_range().end();
    // For each simplex
    for ( ; it_s != it_s_end ; ++it_s)
    {
      Simplex c = *it_s;
      
      int color_simplex = -1;// -1=no color, 0=yellow, 1=red, 2=green, 3=blue
      if (p_simpl_to_color_in_red && 
          std::find(
            p_simpl_to_color_in_red->begin(),
            p_simpl_to_color_in_red->end(),
            c) != p_simpl_to_color_in_red->end())
      {
        color_simplex = 1;
      }
      else if (p_simpl_to_color_in_green && 
          std::find(
            p_simpl_to_color_in_green->begin(),
            p_simpl_to_color_in_green->end(),
            c) != p_simpl_to_color_in_green->end())
      {
        color_simplex = 2;
      }
      else if (p_simpl_to_color_in_blue && 
          std::find(
            p_simpl_to_color_in_blue->begin(),
            p_simpl_to_color_in_blue->end(),
            c) != p_simpl_to_color_in_blue->end())
      {
        color_simplex = 3;
      }

      // Gather the triangles here
      typedef std::vector<Simplex> Triangles;
      Triangles triangles;

      std::size_t num_vertices = c.size();
      // Do not export smaller dimension simplices
      if (num_vertices < m_intrinsic_dimension + 1)
        continue;

      // If m_intrinsic_dimension = 1, each point is output two times,
      // so we need to multiply each index by 2
      // And if only 2 vertices, add a third one (each vertex is duplicated in
      // the file when m_intrinsic dim = 2)
      if (m_intrinsic_dimension == 1)
      {
        std::set<std::size_t> tmp_c;
        std::set<std::size_t>::iterator it = c.begin();
        for ( ; it != c.end() ; ++it)
          tmp_c.insert(*it * 2);
        if (num_vertices == 2)
          tmp_c.insert(*tmp_c.rbegin() + 1);

        c = tmp_c;
      }

      if (num_vertices <= 3)
      {
        triangles.push_back(c);
      }
      else
      {
        // num_vertices >= 4: decompose the simplex in triangles
        std::vector<bool> booleans(num_vertices, false);
        std::fill(booleans.begin() + num_vertices - 3, booleans.end(), true);
        do
        {
          std::set<std::size_t> triangle;
          std::set<std::size_t>::iterator it = c.begin();
          for (int i = 0; it != c.end() ; ++i, ++it)
          {
            if (booleans[i])
              triangle.insert(*it);
          }
          triangles.push_back(triangle);
        }
        while (std::next_permutation(booleans.begin(), booleans.end()));
      }

      // For each cell
      Triangles::const_iterator it_tri = triangles.begin();
      Triangles::const_iterator it_tri_end = triangles.end();
      for ( ; it_tri != it_tri_end ; ++it_tri)
      {
        // Don't export infinite cells
        if (*it_tri->rbegin() == std::numeric_limits<std::size_t>::max())
          continue;

        os << 3 << " ";
        std::set<std::size_t>::const_iterator it_point_idx = it_tri->begin();
        for ( ; it_point_idx != it_tri->end() ; ++it_point_idx)
        {
          os << *it_point_idx << " ";
        }

        if (p_simpl_to_color_in_red || p_simpl_to_color_in_green 
            || p_simpl_to_color_in_blue)
        {
          switch (color_simplex)
          {
            case 0: os << " 255 255 0"; break;
            case 1: os << " 255 0 0"; break;
            case 2: os << " 0 255 0"; break;
            case 3: os << " 0 0 255"; break;
            default: os << " 128 128 128"; break;
          }
        }

        ++num_simplices;
        os << std::endl;
      }
    }

#ifdef CGAL_TC_VERBOSE
    std::cerr << std::endl
      << "=========================================================="
      << std::endl
      << "Export to OFF:\n"
      << "  * Number of vertices: " << m_points.size() << std::endl
      << "  * Total number of simplices: " << num_simplices << std::endl
      << "=========================================================="
      << std::endl;
#endif

    return os;
  }

private:
  const Kernel              m_k;
  const int                 m_intrinsic_dimension;
  const double              m_half_sparsity;
  const double              m_sq_half_sparsity;
  const int                 m_ambient_dim;

  Points                    m_points;
  Weights                   m_weights;
#ifdef CGAL_TC_PERTURB_WEIGHT
  Weights_memory            m_weights_memory;
#endif
#ifdef CGAL_TC_PERTURB_POSITION
  Translations_for_perturb  m_translations;
# if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
  Mutex_for_perturb        *m_p_perturb_mutexes;
# endif
#endif
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
  std::vector<Atomic_wrapper<bool> > m_perturb_tangent_space;
#endif

  Points_ds                 m_points_ds;
  std::vector<bool>         m_are_tangent_spaces_computed;
  TS_container              m_tangent_spaces;
#if defined(CGAL_ALPHA_TC) || defined(CGAL_TC_EXPORT_NORMALS)
  OS_container              m_orth_spaces;
#endif
  Tr_container              m_triangulations; // Contains the triangulations
                                              // and their center vertex
  Stars_container           m_stars;
#ifdef CGAL_LINKED_WITH_TBB
  //std::vector<Tr_mutex>     m_tr_mutexes;
#endif

#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  Points                    m_points_for_tse;
  Points_ds                 m_points_ds_for_tse;
#endif

  mutable CGAL::Random      m_random_generator;

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
