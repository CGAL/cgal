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
#include "CGAL/Tangential_complex/console_color.h"

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/function_objects.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

# include <CGAL/Mesh_3/Profiling_tools.h>

#include <CGAL/IO/Triangulation_off_ostream.h> // CJTODO DEBUG?

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <boost/optional.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/math/special_functions/factorials.hpp>

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

// choose exact integral type for QP solver
// (Gmpzf is not thread-safe)
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
//#define CGAL_QP_NO_ASSERTIONS // CJTODO: NECESSARY? http://doc.cgal.org/latest/QP_solver/group__PkgQPSolverFunctions.html#ga1fefbd0436aca0e281f88e8e6cd8eb74

//#define CGAL_TC_EXPORT_NORMALS // Only for 3D surfaces (k=2, d=3)

//static std::ofstream csv_stream("output/stats.csv"); // CJTODO DEBUG

//CJTODO: debug
//#define CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_SPHERE_2
//#define CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_SPHERE_3
//#define CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_TORUS_D
//#define CGAL_TC_ADD_NOISE_TO_TANGENT_SPACE
//#define CGAL_TC_BETTER_EXPORT_FOR_FLAT_TORUS

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
  typename Kernel_, // ambiant kernel
  typename DimensionTag, // intrinsic dimension
  typename Concurrency_tag = CGAL::Parallel_tag,
#ifdef CGAL_ALPHA_TC
  // For the alpha-TC, the dimension of the RT is variable
  // => we need to force 
  typename Tr = Regular_triangulation
  <
    Epick_d<CGAL::Dynamic_dimension_tag>,
    Triangulation_data_structure
    <
      typename Regular_triangulation_euclidean_traits<
        Epick_d<CGAL::Dynamic_dimension_tag> >::Dimension,
      Triangulation_vertex<Regular_triangulation_euclidean_traits<
        Epick_d<CGAL::Dynamic_dimension_tag> >, Vertex_data >,
      Triangulation_full_cell<Regular_triangulation_euclidean_traits<
        Epick_d<CGAL::Dynamic_dimension_tag> > >
    >
  >
#else
  typename Tr = Regular_triangulation
  <
    Epick_d<DimensionTag>,
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
#endif
>
class Tangential_complex
{
  typedef Kernel_                                     K;
  typedef typename K::FT                              FT;
  typedef typename K::Point_d                         Point;
  typedef typename K::Weighted_point_d                Weighted_point;
  typedef typename K::Vector_d                        Vector;

  typedef Tr                                          Triangulation;
  typedef typename Triangulation::Geom_traits         Tr_traits;
  typedef typename Triangulation::Weighted_point      Tr_point;
  typedef typename Triangulation::Bare_point          Tr_bare_point;
  typedef typename Triangulation::Vertex_handle       Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle    Tr_full_cell_handle;
  typedef typename Tr_traits::Vector_d                Tr_vector;

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

  typedef Point_cloud_data_structure<K, Points>       Points_ds;
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
  typedef Basis<K>                                      Tangent_space_basis;
  typedef Basis<K>                                      Orthogonal_space_basis;
  typedef std::vector<Tangent_space_basis>              TS_container;
  typedef std::vector<Orthogonal_space_basis>           OS_container;

private:
  typedef std::vector<Tr_and_VH>                        Tr_container;
  typedef std::vector<Vector>                           Vectors;

  // An Incident_simplex is the list of the vertex indices
  // except the center vertex
  typedef std::set<std::size_t>                         Incident_simplex;
  typedef std::set<std::size_t>                         Indexed_simplex;
  typedef std::vector<Incident_simplex>                 Star;
  typedef std::vector<Star>                             Stars_container;

  // For the priority queues of solve_inconsistencies_using_alpha_TC
  struct Simplex_and_alpha
  {
    Simplex_and_alpha() {}
    Simplex_and_alpha(
      std::size_t center_point_index,
      Incident_simplex const& simplex, // NOT including "center_point_index"
      FT squared_alpha,
      Vector const& thickening_vector)
    : m_center_point_index(center_point_index),
      m_simplex(simplex),
      m_squared_alpha(squared_alpha),
      m_thickening_vector(thickening_vector)
    {}

    // For the priority queue
    bool operator>(Simplex_and_alpha const& other) const
    {
      return m_squared_alpha > other.m_squared_alpha;
    }

    std::size_t           m_center_point_index; // P
    Incident_simplex      m_simplex; // Missing simplex (does NOT includes P)
    FT                    m_squared_alpha;
    Vector                m_thickening_vector; // (P, Cq)
  };

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
  struct First_of_pair
  {
    template<typename> struct result;

    template <typename F, typename Pair>
    struct result<F(Pair)>
    {
      typedef typename boost::remove_reference<Pair>::type::first_type const& type;
    };

    template <typename Pair>
    typename Pair::first_type const& operator()(Pair const& pair) const
    {
      return pair.first;
    }
  };

public:
  typedef Tangential_complex_::Simplicial_complex     Simplicial_complex;

  /// Constructor for a range of points
  template <typename InputIterator>
  Tangential_complex(InputIterator first, InputIterator last,
                     double sparsity, int intrinsic_dimension,
#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
                     InputIterator first_for_tse, InputIterator last_for_tse,
#endif
                     const K &k = K()
                     )
  : m_k(k),
    m_intrinsic_dim(intrinsic_dimension),
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
#ifdef CGAL_TC_EXPORT_NORMALS
    , m_orth_spaces(m_points.size(), Orthogonal_space_basis())
#endif
#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
    , m_points_for_tse(first_for_tse, last_for_tse)
    , m_points_ds_for_tse(m_points_for_tse)
#endif
  {
    if (sparsity <= 0.)
      std::cerr << "!Warning! Sparsity should be > 0\n";
  }

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
    return m_intrinsic_dim;
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
#ifdef CGAL_TC_EXPORT_NORMALS
                        , const OS_container& orthogonal_spaces
#endif
                         )
  {
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    std::cerr << "Cannot use CGAL_TC_PERTURB_TANGENT_SPACE and set "
              << " tangent spaces manually at the same time\n";
    std::exit(EXIT_FAILURE);
#endif
#ifdef CGAL_TC_EXPORT_NORMALS
    CGAL_assertion(m_points.size() == tangent_spaces.size()
                   && m_points.size() == orthogonal_spaces.size());
#else
    CGAL_assertion(m_points.size() == tangent_spaces.size());
#endif
    m_tangent_spaces = tangent_spaces;
#ifdef CGAL_TC_EXPORT_NORMALS
    m_orth_spaces = orthogonal_spaces;
#endif
    for(std::size_t i=0; i<m_points.size(); ++i)
      m_are_tangent_spaces_computed[i] = true;
  }

  void compute_tangential_complex()
  {
    // CJTODO TEMP
    /*{
    INS_range ins_range = m_points_ds.query_incremental_ANN(
      m_k.construct_point_d_object()(8, ORIGIN));
    int c = 0;
    for (INS_iterator nn_it = ins_range.begin() ;
      nn_it != ins_range.end() && c < 10 ;
      ++nn_it, ++c)
    {
      std::size_t neighbor_point_idx = nn_it->first;
      FT sqdist = nn_it->second;
      std::cerr << neighbor_point_idx << ":" << CGAL::sqrt(sqdist) << "\n";
    }
    }
    {
    INS_range ins_range = m_points_ds.query_incremental_ANN(
      m_points[234]);
    int c = 0;
    for (INS_iterator nn_it = ins_range.begin() ;
      nn_it != ins_range.end() && c < 10 ;
      ++nn_it, ++c)
    {
      std::size_t neighbor_point_idx = nn_it->first;
      FT sqdist = nn_it->second;
      std::cerr << neighbor_point_idx << ":" << CGAL::sqrt(sqdist) << "\n";
    }
    }*/

#ifdef CGAL_TC_PERFORM_EXTRA_CHECKS
    std::cerr << red << "WARNING: CGAL_TC_PERFORM_EXTRA_CHECKS is defined. "
      << "Computation might be slower than usual.\n" << white;
#endif
#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
    std::cerr << red << "WARNING: CGAL_TC_ALVAREZ_SURFACE_WINDOW is defined ("
      << CGAL_TC_ALVAREZ_SURFACE_WINDOW << ").\n" << white;
#endif

#if defined(CGAL_TC_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    Wall_clock_timer t;
#endif

    // We need to do that because we don't want the container to copy the
    // already-computed triangulations (while resizing) since it would
    // invalidate the vertex handles stored beside the triangulations
    m_triangulations.resize(m_points.size());
    m_stars.resize(m_points.size());
    m_squared_star_spheres_radii_incl_margin.resize(m_points.size(), FT(-1));
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
        compute_tangent_triangulation(i);
    }

#if defined(CGAL_TC_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    std::cerr << "Tangential complex computed in " << t.elapsed()
              << " seconds.\n";
#endif
  }

  void estimate_intrinsic_dimension() const
  {
    // Kernel functors
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

    std::vector<FT> sum_eigen_values(m_ambient_dim, FT(0));
    std::size_t num_points_for_pca = static_cast<std::size_t>(
      std::pow(BASE_VALUE_FOR_PCA, m_intrinsic_dim));

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
#if defined(CGAL_TC_VERBOSE) || defined(CGAL_TC_PROFILING)
    std::cerr << yellow << "\nRefreshing TC... " << white;
#endif

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
        compute_tangent_triangulation(i);
    }

#ifdef CGAL_TC_PROFILING
    std::cerr << yellow << "done in " << t.elapsed()
      << " seconds.\n" << white;
#elif defined(CGAL_TC_VERBOSE)
    std::cerr << yellow << "done.\n" << white;
#endif
  }

  // If the list of perturbed points is provided, it is much faster
  template <typename Point_indices_range>
  void refresh_tangential_complex(
    Point_indices_range const& perturbed_points_indices)
  {
#if defined(CGAL_TC_VERBOSE) || defined(CGAL_TC_PROFILING)
    std::cerr << yellow << "\nRefreshing TC... " << white;
#endif

#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif

    // ANN tree containing only the perturbed points
    Points_ds updated_pts_ds(m_points, perturbed_points_indices);

#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()),
        Refresh_tangent_triangulation(*this, updated_pts_ds)
        );
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      for (std::size_t i = 0 ; i < m_points.size() ; ++i)
        refresh_tangent_triangulation(i, updated_pts_ds);
    }

#ifdef CGAL_TC_PROFILING
    std::cerr << yellow << "done in " << t.elapsed()
      << " seconds.\n" << white;
#elif defined(CGAL_TC_VERBOSE)
    std::cerr << yellow << "done.\n" << white;
#endif
  }

  // time_limit in seconds: 0 = no fix to do, < 0 = no time limit
  Fix_inconsistencies_status fix_inconsistencies_using_perturbation(
    unsigned int &num_steps,
    std::size_t &initial_num_inconsistent_stars,
    std::size_t &best_num_inconsistent_stars,
    std::size_t &final_num_inconsistent_stars,
    double time_limit = -1.)
  {
    if (time_limit == 0.)
      return TIME_LIMIT_REACHED;

    Wall_clock_timer t;

#ifdef CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
    CGAL::cpp11::tuple<std::size_t, std::size_t, std::size_t> stats_before =
      number_of_inconsistent_simplices(false);

    if (CGAL::cpp11::get<1>(stats_before) == 0)
    {
# ifdef CGAL_TC_VERBOSE
      std::cerr << "Nothing to fix.\n";
# endif
      return TC_FIXED;
    }
#endif // CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

    bool done = false;
    best_num_inconsistent_stars = m_triangulations.size();
    num_steps = 0;
    while (!done)
    {
#ifdef CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
      std::cerr
        << "\nBefore fix step:\n"
        << "  * Total number of simplices in stars (incl. duplicates): "
        << CGAL::cpp11::get<0>(stats_before) << "\n"
        << "  * Num inconsistent simplices in stars (incl. duplicates): "
        << red << CGAL::cpp11::get<1>(stats_before) << white << " (" 
        << 100. * CGAL::cpp11::get<1>(stats_before) / CGAL::cpp11::get<0>(stats_before) << "%)\n"
        << "  * Number of stars containing inconsistent simplices: "
        << red << CGAL::cpp11::get<2>(stats_before) << white << " ("
        << 100. * CGAL::cpp11::get<2>(stats_before) / m_points.size() << "%)\n";
#endif

#if defined(CGAL_TC_VERBOSE) || defined(CGAL_TC_PROFILING)
      std::cerr << yellow 
        << "\nAttempt to fix inconsistencies using perturbations - step #"
        << num_steps + 1 << "... " << white;
#endif

      std::size_t num_inconsistent_stars = 0;
      std::vector<std::size_t> updated_points;

#ifdef CGAL_TC_PROFILING
      Wall_clock_timer t_fix_step;
#endif

      // Parallel
#if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
      if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
      {
        tbb::combinable<std::size_t> num_inconsistencies;
        tbb::combinable<std::vector<std::size_t> > tls_updated_points;
        tbb::parallel_for(
          tbb::blocked_range<size_t>(0, m_triangulations.size()),
          Try_to_solve_inconsistencies_in_a_local_triangulation(
          *this, num_inconsistencies, tls_updated_points)
        );
        num_inconsistent_stars =
          num_inconsistencies.combine(std::plus<std::size_t>());
        updated_points = tls_updated_points.combine(
          [](std::vector<std::size_t> const& x, std::vector<std::size_t> const& y) { // CJTODO: C++11
            std::vector<std::size_t> res;
            res.reserve(x.size() + y.size());
            res.insert(res.end(), x.begin(), x.end());
            res.insert(res.end(), y.begin(), y.end());
            return res;
          });
      }
      // Sequential
      else
#endif // CGAL_LINKED_WITH_TBB
      {
        for (std::size_t i = 0 ; i < m_triangulations.size() ; ++i)
        {
          num_inconsistent_stars +=
            try_to_solve_inconsistencies_in_a_local_triangulation(
              i, std::back_inserter(updated_points));
        }
      }

      double fix_step_time = t_fix_step.elapsed();

#if defined(CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES) || defined(CGAL_TC_VERBOSE)
      std::cerr
        << "\nEncountered during fix:\n"
        << "  * Num stars containing inconsistent simplices: "
        << red << num_inconsistent_stars << white
        << " (" << 100. * num_inconsistent_stars / m_points.size() << "%)\n";
#endif

#ifdef CGAL_TC_PROFILING
      std::cerr << yellow << "done in " << fix_step_time
        << " seconds.\n" << white;
#elif defined(CGAL_TC_VERBOSE)
      std::cerr << yellow << "done.\n" << white;
#endif

#ifdef CGAL_TC_GLOBAL_REFRESH
      if (num_inconsistent_stars > 0)
        refresh_tangential_complex(updated_points);

# ifdef CGAL_TC_PERFORM_EXTRA_CHECKS
      // Confirm that all stars were actually refreshed
      std::size_t num_inc_1 = 
        CGAL::cpp11::get<1>(number_of_inconsistent_simplices(false));
      refresh_tangential_complex();
      std::size_t num_inc_2 = 
        CGAL::cpp11::get<1>(number_of_inconsistent_simplices(false));
      if (num_inc_1 != num_inc_2)
        std::cerr << red << "REFRESHMENT CHECK: FAILED. ("
          << num_inc_1 << " vs " << num_inc_2 << ")\n" << white;
      else
        std::cerr << green << "REFRESHMENT CHECK: PASSED.\n" << white;
# endif
#endif

#ifdef CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
      CGAL::cpp11::tuple<std::size_t, std::size_t, std::size_t> stats_after =
        number_of_inconsistent_simplices(false);

      std::cerr
        << "\nAfter fix:\n"
        << "  * Total number of simplices in stars (incl. duplicates): "
        << CGAL::cpp11::get<0>(stats_after) << "\n"
        << "  * Num inconsistent simplices in stars (incl. duplicates): "
        << red << CGAL::cpp11::get<1>(stats_after) << white << " (" 
        << 100. * CGAL::cpp11::get<1>(stats_after) / CGAL::cpp11::get<0>(stats_after) << "%)\n"
        << "  * Number of stars containing inconsistent simplices: "
        << red << CGAL::cpp11::get<2>(stats_after) << white << " ("
        << 100. * CGAL::cpp11::get<2>(stats_after) / m_points.size() << "%)\n";

      stats_before = stats_after;
#endif

      if (num_steps == 0)
        initial_num_inconsistent_stars = num_inconsistent_stars;

      if (num_inconsistent_stars < best_num_inconsistent_stars)
        best_num_inconsistent_stars = num_inconsistent_stars;

      final_num_inconsistent_stars = num_inconsistent_stars;

      done = (num_inconsistent_stars == 0);
      if (!done)
      {
        ++num_steps;
        if (time_limit > 0. && t.elapsed() > time_limit)
        {
#ifdef CGAL_TC_VERBOSE
          std::cerr << red << "Time limit reached.\n" << white;
#endif
          return TIME_LIMIT_REACHED;
        }
      }
    }

#ifdef CGAL_TC_VERBOSE
    std::cerr << green << "Fixed!\n" << white;
#endif
    return TC_FIXED;
  }

  // Return a tuple
  // <num_simplices, num_inconsistent_simplices, num_inconsistent_stars>
  CGAL::cpp11::tuple<std::size_t, std::size_t, std::size_t> 
  number_of_inconsistent_simplices(
#ifdef CGAL_TC_VERBOSE
    bool verbose = true
#else
    bool verbose = false
#endif
    ) const
  {
    std::size_t num_simplices = 0;
    std::size_t num_inconsistent_simplices = 0;
    std::size_t num_inconsistent_stars = 0;
    // For each triangulation
    for (std::size_t idx = 0 ; idx < m_points.size() ; ++idx)
    {
      bool is_star_inconsistent = false;

      // For each cell
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        // Don't check infinite cells
        if (is_infinite(*it_inc_simplex))
          continue;

        Indexed_simplex c = *it_inc_simplex;
        c.insert(idx); // Add the missing index

        if (!is_simplex_consistent(c))
        {
          ++num_inconsistent_simplices;
          is_star_inconsistent = true;
        }

        ++num_simplices;
      }
      num_inconsistent_stars += is_star_inconsistent;
    }

    if (verbose)
    {
      std::cerr
        << "\n==========================================================\n"
        << "Inconsistencies:\n"
        << "  * Total number of simplices in stars (incl. duplicates): "
        << num_simplices << "\n"
        << "  * Number of inconsistent simplices in stars (incl. duplicates): "
        << num_inconsistent_simplices << " (" 
        << 100. * num_inconsistent_simplices / num_simplices << "%)\n"
        << "  * Number of stars containing inconsistent simplices: "
        << num_inconsistent_stars << " ("
        << 100. * num_inconsistent_stars / m_points.size() << "%)\n"
        << "==========================================================\n";
    }

    return std::make_tuple(
      num_simplices, num_inconsistent_simplices, num_inconsistent_stars);
  }

  // First clears the complex then exports the TC into it
  // Return the max dimension of the simplices
  // check_lower_and_higher_dim_simplices : 0 (false), 1 (true), 2 (auto)
  //   If the check is enabled, the function:
  //   - won't insert the simplex if it is already in a higher dim simplex
  //   - will erase any lower-dim simplices that are faces of the new simplex
  //   "auto" (= 2) will enable the check as a soon as it encounters a 
  //   simplex whose dimension is different from the previous ones.
  //   N.B.: The check is quite expensive.
  int export_TC(Simplicial_complex &complex,
    bool export_infinite_simplices = false, 
    int check_lower_and_higher_dim_simplices = 2,
    std::set<Indexed_simplex > *p_inconsistent_simplices = NULL) const
  {
#if defined(CGAL_TC_VERBOSE) || defined(CGAL_TC_PROFILING)
    std::cerr << yellow
      << "\nExporting the TC as a Simplicial_complex... " << white;
#endif
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif

    int max_dim = -1;
    complex.clear();

    // For each triangulation
    for (std::size_t idx = 0 ; idx < m_points.size() ; ++idx)
    {
      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        Indexed_simplex c = *it_inc_simplex;

        // Don't export infinite cells
        if (!export_infinite_simplices && is_infinite(c))
          continue;

        // Unusual simplex dim?
        if (check_lower_and_higher_dim_simplices == 2
          && max_dim != -1
          && static_cast<int>(c.size()) != max_dim)
        {
          // Let's activate the check
          std::cerr << red << 
            "Info: check_lower_and_higher_dim_simplices ACTIVATED. "
            "Export might be take some time...\n" << white;
          check_lower_and_higher_dim_simplices = 1;
        }

        if (static_cast<int>(c.size()) > max_dim)
          max_dim = static_cast<int>(c.size());
        // Add the missing center vertex
        c.insert(idx);

#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
        if (is_one_of_the_coord_far_from_origin(c, CGAL_TC_ALVAREZ_SURFACE_WINDOW, 2))
          continue;
#endif

        // Try to insert the simplex
        bool added = 
          complex.add_simplex(c, check_lower_and_higher_dim_simplices == 1);

        // Inconsistent?
        if (p_inconsistent_simplices && added && !is_simplex_consistent(c))
        {
          p_inconsistent_simplices->insert(c);
        }
      }
    }

#ifdef CGAL_TC_PROFILING
    std::cerr << yellow << "done in " << t.elapsed()
      << " seconds.\n" << white;
#elif defined(CGAL_TC_VERBOSE)
    std::cerr << yellow << "done.\n" << white;
#endif

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
  }

  // Returns true if some inconsistencies were found
  bool check_and_solve_inconsistencies_by_filtering_simplices_out()
  {
    bool inconsistencies_found = false;

    // CJTODO: parallel_for???
    for (std::size_t idx = 0 ; idx < m_triangulations.size() ; ++idx)
    {
      if (filter_inconsistent_simplices_in_a_local_triangulation(idx))
        inconsistencies_found = true;
    }

    return inconsistencies_found;
  }


#ifdef CGAL_ALPHA_TC
private:
  
  // Look in the star of point "i" for inconsistent simplices, compute
  // an approximation of alpha for each one and push it into the priority
  // queues.
  // Returns the number of inconsistent simplices found
  template <typename PQueues>
  std::size_t fill_pqueues_for_alpha_tc(std::size_t i, PQueues &pqueues)
  {
    // Kernel/traits functors
    typename K::Difference_of_points_d k_diff_points =
      m_k.difference_of_points_d_object();
    typename K::Squared_length_d k_sqlen =
      m_k.squared_length_d_object();
    typename Tr_traits::Construct_weighted_point_d constr_wp =
      m_triangulations[0].tr().geom_traits().construct_weighted_point_d_object();

    std::size_t num_inconsistent_simplices = 0;

    // For each cell
    Star::const_iterator it_inc_simplex = m_stars[i].begin();
    Star::const_iterator it_inc_simplex_end = m_stars[i].end();
    for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
    {
      Incident_simplex const& s = *it_inc_simplex;

      // Don't check infinite cells
      if (is_infinite(s))
        continue;

      int simplex_dim = static_cast<int>(s.size());

      // P: points whose star does not contain "s"
      std::vector<std::size_t> P;
      is_simplex_consistent(i, s, std::back_inserter(P), true);

      if (!P.empty())
      {
        ++num_inconsistent_simplices;

        Triangulation const& q_tr = m_triangulations[i].tr();
        Indexed_simplex full_simplex = s;
        full_simplex.insert(i);
        for (std::vector<std::size_t>::const_iterator it_p = P.begin(),
          it_p_end = P.end() ; it_p != it_p_end ; ++it_p)
        {
          // star(p) does not contain "s"
          std::size_t p = *it_p;

          boost::optional<Tr_bare_point> intersection;

          // If we know the intersection between Tq and voronoi(s) is a point
          // I.e. dim(Tq) + dim(V(full_simplex)) = D
          //  <=> dim(Tq) = dim(full_simplex)
          if (q_tr.current_dimension() == s.size())
          {
            // Compute the intersection between aff(Voronoi_cell(s)) and Tq
            // Note that the computation in done in the sub-space defined by Tq
            intersection = 
              compute_aff_of_voronoi_face_and_tangent_subspace_intersection(
                q_tr.current_dimension(),
                project_points_and_compute_weights(
                  full_simplex, m_tangent_spaces[i], q_tr.geom_traits()),
                m_tangent_spaces[i],
                q_tr.geom_traits());

            // CJTODO: replace with an assertion?
            if (!intersection)
            {
              std::cerr << "ERROR fill_pqueues_for_alpha_tc: "
                "aff(Voronoi_cell(s)) and Tq do not intersect.\n";
              continue;
            }
          }
          else
          {
            std::cerr << "***************** Intersection is not a point ********************\n"; // CJTODO DEBUG
            // Compute the intersection between Voronoi_cell(s) and Tq
            // We need Q, i.e. the vertices which are common neighbors of all 
            // vertices of full_simplex in the Delaunay triangulation, but we 
            // don't have, so we use the 5^d nearest neighbors instead
            int num_neighbors_for_Q_approx = std::pow(
              BASE_VALUE_FOR_ALPHA_TC_NEIGHBORHOOD, m_intrinsic_dim);
            const Point center_pt = compute_perturbed_point(p);
            KNS_range ins_range = m_points_ds.query_ANN(
              center_pt, num_neighbors_for_Q_approx);
            
            // CJTODO: optimize that, use a vector!
            std::set<std::size_t> Q(
              boost::make_counting_iterator(std::size_t(0)),
              boost::make_counting_iterator(m_points.size()));
            for (auto i : full_simplex)
              Q.erase(i);

            intersection =
              compute_voronoi_face_and_tangent_subspace_intersection(
              q_tr.current_dimension(),
              project_points_and_compute_weights(
              full_simplex, m_tangent_spaces[i], q_tr.geom_traits()),
              project_points_and_compute_weights(
                Q,
                //boost::counting_range(std::size_t(0), m_points.size() - 1),
                //boost::adaptors::transform(ins_range, First_of_pair()),
                m_tangent_spaces[i], q_tr.geom_traits()),
              m_tangent_spaces[i],
              q_tr.geom_traits());

            // CJTODO: replace with an assertion?
            if (!intersection)
            {
              std::cerr << "ERROR fill_pqueues_for_alpha_tc: "
                "Voronoi_cell(s) and Tq do not intersect.\n";
              continue;
            }
          }

          // The following computations are done in the Euclidian space
          Point inters_global = unproject_point(
            constr_wp(*intersection, 0), m_tangent_spaces[i], 
            q_tr.geom_traits());
          Vector thickening_v = k_diff_points(
            inters_global, compute_perturbed_point(p));
          FT squared_alpha = k_sqlen(thickening_v);

          // We insert full_simplex \ p
          Incident_simplex is = full_simplex;
          is.erase(p);
          pqueues[simplex_dim - m_intrinsic_dim].push(
            Simplex_and_alpha(p, is, squared_alpha, thickening_v));

          // CJTODO DEBUG
          /*std::cerr 
            << "Just inserted the simplex ";
          std::copy(full_simplex.begin(), full_simplex.end(), 
            std::ostream_iterator<std::size_t>(std::cerr, ", ")); 
          std::cerr << "into pqueue (i = " << i << ")\n";*/
        }
      }
    }

    return num_inconsistent_simplices;
  }

public:
  void solve_inconsistencies_using_alpha_TC()
  {
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif

#ifdef CGAL_TC_VERBOSE
    std::cerr << "Fixing inconsistencies using alpha TC...\n";
#endif

    //-------------------------------------------------------------------------
    // 1. Fill priority queues
    //-------------------------------------------------------------------------
    typedef std::priority_queue<Simplex_and_alpha,
      std::vector<Simplex_and_alpha>,
      std::greater<Simplex_and_alpha> > AATC_pq;
    typedef std::vector<AATC_pq> PQueues;
    
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t_pq;
#endif

    // One queue per dimension, from intrinsic dim (index = 0) to 
    // ambiant dim (index = ambiant - intrinsic dim)
    PQueues pqueues;
    pqueues.resize(m_ambient_dim - m_intrinsic_dim + 1);

    std::size_t num_inconsistent_simplices = 0;
    // For each triangulation
    for (std::size_t i = 0 ; i < m_points.size() ; ++i)
      num_inconsistent_simplices += fill_pqueues_for_alpha_tc(i, pqueues);

#ifdef CGAL_TC_VERBOSE
    std::cerr
      << "Num inconsistent simplices found when filling the priority queues: "
      << num_inconsistent_simplices;
# ifdef CGAL_TC_PROFILING
    std::cerr << " (" << t_pq.elapsed() << " s)\n";
# endif
    std::cerr << "\n";
#endif

    //-------------------------------------------------------------------------
    // 2. Thicken tangent spaces to solve inconsistencies
    //-------------------------------------------------------------------------
    
    // While there's elements to treat in the queues
    for(;;)
    {
#ifdef CGAL_TC_PROFILING
      Wall_clock_timer t_one_fix;
#endif

      // Pick the simplex with the lowest dimension and the lowest alpha
      Simplex_and_alpha saa;
      bool found_saa = false;
      for (PQueues::iterator it_pq = pqueues.begin(), it_pq_end = pqueues.end();
        !found_saa && it_pq != it_pq_end ; ++it_pq)
      {
        while (!found_saa && !it_pq->empty())
        {
          saa = it_pq->top();
          it_pq->pop();

          // Check if the simplex is still missing in the star
          if (!is_simplex_in_star(saa.m_center_point_index, saa.m_simplex))
            found_saa = true;
        }
      }

      // If all the queues are empty, we're done!
      if (!found_saa)
        break;

      Tangent_space_basis &tangent_basis = 
        m_tangent_spaces[saa.m_center_point_index];

      // If we're already in the ambiant dim, we just need to thicken the
      // tangent subspace a bit more (see below)
      if (tangent_basis.dimension() < m_ambient_dim)
      {
        // Otherwise, let's thicken the tangent space
        bool vec_added = add_vector_to_orthonormal_basis(
          tangent_basis,
          saa.m_thickening_vector,
          m_k,
          FT(0), /* sqlen_threshold: default value */
          true /* add_to_thickening_vectors */);

        // CJTODO DEBUG
        if (!vec_added)
        {
          std::cerr << "FYI: the thickening vector was not added "
                       "to the basis since it was linearly dependent to it.\n";
        }
      }

      // Update the alpha+/- values
      tangent_basis.update_alpha_values_of_thickening_vectors(
        saa.m_thickening_vector, m_k);

#ifdef CGAL_TC_PROFILING
      Wall_clock_timer t_recomputation;
#endif
      // Recompute triangulation & star
      compute_tangent_triangulation(saa.m_center_point_index);
#ifdef CGAL_TC_PROFILING
      double recomp_timing = t_recomputation.elapsed();
#endif


#ifdef CGAL_TC_PERFORM_EXTRA_CHECKS
      if (!is_simplex_in_star(saa.m_center_point_index, saa.m_simplex))
      {
        std::cerr 
          << "FAILED in solve_inconsistencies_using_alpha_TC(): "
          << "simplex " << saa.m_center_point_index << ", ";
        std::copy(saa.m_simplex.begin(), saa.m_simplex.end(), 
          std::ostream_iterator<std::size_t>(std::cerr, ", ")); 
        std::cerr << " not added in star #"
          << saa.m_center_point_index
          << " (basis dim = " << tangent_basis.dimension()
# ifdef CGAL_TC_PROFILING
          << " - " << t_one_fix.elapsed() << " s [recomputation = "
          << recomp_timing << " s]"
# endif
          << ")\n";

        Indexed_simplex full_s = saa.m_simplex;
        full_s.insert(saa.m_center_point_index);

        // CJTODO DEBUG
        bool is_this_simplex_somewhere = false;
        for(auto ii : saa.m_simplex) // CJTODO C++11
        {
          Indexed_simplex z = full_s;
          z.erase(ii);
          if (is_simplex_in_star(ii, z))
          {
            is_this_simplex_somewhere = true;
            std::cerr << "The simplex is in star #" << ii << "\n";
            break;
          }
        }
        if (!is_this_simplex_somewhere)
          std::cerr << "WOW The simplex is NOWHERE!\n";

        // CJTODO DEBUG
        if (m_ambient_dim <= 3)
        {
          if (is_simplex_in_the_ambient_delaunay(full_s))
            std::cerr << "The simplex is in the ambiant Delaunay.\n";
          else
            std::cerr << "The simplex is NOT in the ambiant Delaunay.\n";

          std::cerr << "Checking simplices of the star #" 
            << saa.m_center_point_index << "\n";
          Star const& star = m_stars[saa.m_center_point_index];
          for (Star::const_iterator is = star.begin(), is_end = star.end() ;
            is != is_end ; ++is)
          { 
            if (is_simplex_in_the_ambient_delaunay(*is))
              std::cerr << "The simplex is in the ambiant Delaunay.\n";
            else
            {
              std::cerr << "The simplex is NOT in the ambiant Delaunay.\n";
              for(auto ii : *is) // CJTODO C++11
                perturb(ii);
            }
          }
        }
         
        std::cerr << "Perturbing the points...\n";
        perturb(saa.m_center_point_index);
        for(auto ii : saa.m_simplex) // CJTODO C++11
          perturb(ii);
        refresh_tangential_complex();
        pqueues.clear();
        pqueues.resize(m_ambient_dim - m_intrinsic_dim + 1);
        
        std::size_t num_inconsistent_simplices = 0;
        // For each triangulation
        for (std::size_t i = 0 ; i < m_points.size() ; ++i)
          num_inconsistent_simplices += fill_pqueues_for_alpha_tc(i, pqueues);

#ifdef CGAL_TC_VERBOSE
        std::cerr
          << "Num inconsistent simplices found when filling the priority queues: "
          << num_inconsistent_simplices << "\n";
#endif
      }
      // CJTODO DEBUG
      else
      {
        std::cerr << "SUCCESS: "
          << saa.m_center_point_index << ", ";
        std::copy(saa.m_simplex.begin(), saa.m_simplex.end(), 
          std::ostream_iterator<std::size_t>(std::cerr, ", ")); 
        std::cerr << " added in star #" 
          << saa.m_center_point_index 
          << " (basis dim = " << tangent_basis.dimension()
# ifdef CGAL_TC_PROFILING
          << " - " << t_one_fix.elapsed() << " s [recomputation = "
          << recomp_timing << " s]"
# endif
          << ")\n";
        //check_if_all_simplices_are_in_the_ambient_delaunay();
      }
#endif
      
      // It's not a problem if entries are duplicated in the pqueues
      // since there's a check when we pop elements
      fill_pqueues_for_alpha_tc(saa.m_center_point_index, pqueues);
    }

#ifdef CGAL_TC_PROFILING
    std::cerr << "Tangential complex fixed in " << t.elapsed()
              << " seconds.\n";
#endif
  }
#endif // CGAL_ALPHA_TC
  
  template<typename ProjectionFunctor = CGAL::Identity<Point> >
  std::ostream &export_to_off(
    const Simplicial_complex &complex, std::ostream & os,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_red = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_green = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_blue = NULL,
    ProjectionFunctor const& point_projection = ProjectionFunctor())
    const
  {
    return export_to_off(
      os, false, p_simpl_to_color_in_red, p_simpl_to_color_in_green, 
      p_simpl_to_color_in_blue, &complex, point_projection);
  }
  
  template<typename ProjectionFunctor = CGAL::Identity<Point> >
  std::ostream &export_to_off(
    std::ostream & os, bool color_inconsistencies = false,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_red = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_green = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_blue = NULL,
    const Simplicial_complex *p_complex = NULL,
    ProjectionFunctor const& point_projection = ProjectionFunctor()) const
  {
    if (m_points.empty())
      return os;

    if (m_ambient_dim < 2)
    {
      std::cerr << "Error: export_to_off => ambient dimension should be >= 2.\n";
      os << "Error: export_to_off => ambient dimension should be >= 2.\n";
      return os;
    }
    if (m_ambient_dim > 3)
    {
      std::cerr << "Warning: export_to_off => ambient dimension should be "
                   "<= 3. Only the first 3 coordinates will be exported.\n";
    }

    if (m_intrinsic_dim < 1 || m_intrinsic_dim > 3)
    {
      std::cerr << "Error: export_to_off => intrinsic dimension should be "
                   "between 1 and 3.\n";
      os << "Error: export_to_off => intrinsic dimension should be "
            "between 1 and 3.\n";
      return os;
    }

    std::stringstream output;
    std::size_t num_simplices, num_vertices;
    export_vertices_to_off(output, num_vertices, false, point_projection);
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
    // For each triangulation
    for (std::size_t idx = 0 ; idx < m_points.size() ; ++idx)
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
        if (is_infinite(*it_inc_simplex))
          continue;

        Indexed_simplex c = *it_inc_simplex;
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
  
  // Expensive!
  bool is_simplex_in_the_ambient_delaunay(
    Indexed_simplex const& s) const
  {
    //-------------------------------------------------------------------------
    // Build the ambient Delaunay triangulation
    // Then save its simplices into "amb_dt_simplices"
    //-------------------------------------------------------------------------

    typedef Regular_triangulation_euclidean_traits<K>         RT_Traits;
    typedef Regular_triangulation<
      RT_Traits,
      Triangulation_data_structure<
        typename RT_Traits::Dimension,
        Triangulation_vertex<RT_Traits, Vertex_data>
      > >                                                     RT;
    typedef typename RT::Vertex_handle                        RT_VH;
    typedef typename RT::Finite_full_cell_const_iterator      FFC_it;

    RT ambient_dt(m_ambient_dim);
    for (std::size_t i=0; i<m_points.size(); ++i)
    {
      const Weighted_point wp = compute_perturbed_weighted_point(i);
      RT_VH vh = ambient_dt.insert(wp);
      vh->data() = i;
    }

    for (FFC_it cit = ambient_dt.finite_full_cells_begin() ;
         cit != ambient_dt.finite_full_cells_end() ; ++cit )
    {
      Indexed_simplex simplex;
      for (int i = 0 ; i < m_ambient_dim + 1 ; ++i)
        simplex.insert(cit->vertex(i)->data());

      if (std::includes(simplex.begin(), simplex.end(), 
                        s.begin(), s.end()))
        return true;
    }

    return false;
  }

  bool check_if_all_simplices_are_in_the_ambient_delaunay(
    const Simplicial_complex *p_complex = NULL,
    bool check_for_any_dimension_simplices = true,
    std::set<Indexed_simplex > * incorrect_simplices = NULL) const
  {
    typedef Simplicial_complex::Simplex                     Simplex;
    typedef Simplicial_complex::Simplex_set                 Simplex_set;

    if (m_points.empty())
      return true;

    typedef Regular_triangulation_euclidean_traits<K>       RT_Traits;
    typedef Regular_triangulation<
      RT_Traits,
      Triangulation_data_structure<
        typename RT_Traits::Dimension,
        Triangulation_vertex<RT_Traits, Vertex_data>
      > >                                                   RT;
    typedef typename RT::Vertex_handle                      RT_VH;
    typedef typename RT::Finite_full_cell_const_iterator    FFC_it;

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
        (check_for_any_dimension_simplices ? 1 : m_intrinsic_dim);
      int highest_dim =
        (check_for_any_dimension_simplices ? m_ambient_dim : m_intrinsic_dim);

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

    Simplex_set const *p_simplices;

    std::size_t num_infinite_cells = 0;
    Simplex_set stars_simplices;
    if (!p_complex)
    {
      Stars_container::const_iterator it_star = m_stars.begin();
      Stars_container::const_iterator it_star_end = m_stars.end();
      // For each star: get the finite simplices
      for ( ; it_star != it_star_end ; ++it_star)
      {
        for (Star::const_iterator it_s = it_star->begin(), 
          it_s_end = it_star->end() ; it_s != it_s_end ; ++it_s)
        {
          if (!is_infinite(*it_s))
            stars_simplices.insert(*it_s);
        }
      }
      /*typename Tr_container::const_iterator it_tr = m_triangulations.begin();
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
            ++num_infinite_cells;
            continue;
          }
          Simplex simplex;
          for (int i = 0 ; i < tr.current_dimension() + 1 ; ++i)
            simplex.insert((*it_c)->vertex(i)->data());

          stars_simplices.insert(simplex);
        }
      }*/



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

#ifdef CGAL_TC_VERBOSE
    std::cerr
      << (incorrect_simplices->empty() ? "OK " : "ERROR ")
      << "check_if_all_simplices_are_in_the_ambient_delaunay:\n"
      << "  Number of simplices in ambient RT: " << amb_dt_simplices.size()
      << "\n"
      << "  Number of unique simplices in TC stars: " << p_simplices->size()
      << "\n"
      << "  Number of infinite full cells in TC stars: " << num_infinite_cells
      << "\n"
      << "  Number of wrong simplices: " << incorrect_simplices->size()
      << "\n";
#endif
    return incorrect_simplices->empty();
  }

private:

  class Compare_distance_to_ref_point
  {
  public:
    Compare_distance_to_ref_point(Point const& ref, K const& k)
      : m_ref(ref), m_k(k) {}

    bool operator()(Point const& p1, Point const& p2)
    {
      typename K::Squared_distance_d sqdist =
        m_k.squared_distance_d_object();
      return sqdist(p1, m_ref) < sqdist(p2, m_ref);
    }

  private:
    Point const& m_ref;
    K const& m_k;
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
        m_tc.compute_tangent_triangulation(i);
    }
  };

  // Functor for resfresh_tangential_complex function
  class Refresh_tangent_triangulation
  {
    Tangential_complex & m_tc;
    Points_ds const& m_updated_pts_ds;

  public:
    // Constructor
    Refresh_tangent_triangulation(
      Tangential_complex &tc, Points_ds const& updated_pts_ds)
      : m_tc(tc), m_updated_pts_ds(updated_pts_ds)
    { }

    // Constructor
    Refresh_tangent_triangulation(const Refresh_tangent_triangulation &ctt)
      : m_tc(ctt.m_tc), m_updated_pts_ds(ctt.m_updated_pts_ds)
    { }

    // operator()
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
      for (size_t i = r.begin() ; i != r.end() ; ++i)
        m_tc.refresh_tangent_triangulation(i, m_updated_pts_ds);
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  bool is_infinite(Indexed_simplex const& s) const
  {
    return *s.rbegin() == std::numeric_limits<std::size_t>::max();
  }
  
  bool is_one_of_the_coord_far_from_origin(
    Point const& p, FT limit, int only_test_the_first_n_coords = -1) const
  {
    typename K::Construct_cartesian_const_iterator_d ccci =
      m_k.construct_cartesian_const_iterator_d_object();
    int c = 1;
    for (auto it = ccci(p) ; it != ccci(p, 0) ; ++it, ++c) // CJTODO: C++11
    {
      if (*it > limit || *it < -limit)
        return true;

      if (only_test_the_first_n_coords > 0 && c == only_test_the_first_n_coords)
        break;
    }
    return false;
  }

  bool is_one_of_the_coord_far_from_origin(
    Indexed_simplex const& s, FT limit, int only_test_the_first_n_coords = -1) const
  {
    for (Indexed_simplex::const_iterator it_index = s.begin();
      it_index != s.end() ; ++it_index)
    {
      // Infinite vertex? Much too far!
      if (*it_index == std::numeric_limits<std::size_t>::max())
        return true;

      if (is_one_of_the_coord_far_from_origin(
        compute_perturbed_point(*it_index), limit, only_test_the_first_n_coords))
        return true;
    }
    return false;
  }

  // Output: "triangulation" is a Regular Triangulation containing at least the
  // star of "center_pt"
  // Returns the handle of the center vertex
  Tr_vertex_handle compute_star(
    std::size_t i, const Point &center_pt, const Tangent_space_basis &tsb,
    Triangulation &triangulation, bool verbose = false)
  {
    int tangent_space_dim = tsb.dimension();
    const Tr_traits &local_tr_traits = triangulation.geom_traits();
    Tr_vertex_handle center_vertex;

    // Kernel functor & objects
    typename K::Squared_distance_d k_sqdist = m_k.squared_distance_d_object();

    // Triangulation's traits functor & objects
    typename Tr_traits::Point_weight_d point_weight =
      local_tr_traits.point_weight_d_object();
    typename Tr_traits::Power_center_d power_center =
      local_tr_traits.power_center_d_object();

    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p)
    //***************************************************

    // Insert p
    Tr_point proj_wp;
    if (i == tsb.origin())
    {
      // Insert {(0, 0, 0...), m_weights[i]}
      proj_wp = local_tr_traits.construct_weighted_point_d_object()(
        local_tr_traits.construct_point_d_object()(tangent_space_dim, ORIGIN),
        m_weights[i]);
    }
    else
    {
      const Weighted_point& wp = compute_perturbed_weighted_point(i);
      proj_wp = project_point_and_compute_weight(wp, tsb, local_tr_traits);
    }

    center_vertex = triangulation.insert(proj_wp);
    center_vertex->data() = i;
    if (verbose)
      std::cerr << "* Inserted point #" << i << "\n";

#ifdef CGAL_TC_VERY_VERBOSE
    std::size_t num_attempts_to_insert_points = 1;
    std::size_t num_inserted_points = 1;
#endif
    //const int NUM_NEIGHBORS = 150;
    //KNS_range ins_range = m_points_ds.query_ANN(center_pt, NUM_NEIGHBORS);
    INS_range ins_range = m_points_ds.query_incremental_ANN(center_pt);

    // While building the local triangulation, we keep the radius
    // of the sphere "star sphere" centered at "center_vertex"
    // and which contains all the
    // circumspheres of the star of "center_vertex"
    boost::optional<FT> squared_star_sphere_radius_plus_margin;

#ifdef CGAL_ALPHA_TC  
    /*FT max_absolute_alpha = tsb.max_absolute_alpha();
    // "2*m_half_sparsity" because both points can be perturbed
    FT max_sqdist_to_tangent_space = (max_absolute_alpha == FT(0) ?
    FT(0) : CGAL::square(2*max_absolute_alpha + 2*m_half_sparsity));*/
    std::size_t number_of_attempts_to_insert_points = 0;
    const std::size_t MAX_NUM_INSERTED_POINTS =
      tsb.num_thickening_vectors() > 0 ?
      static_cast<std::size_t>(std::pow(BASE_VALUE_FOR_ALPHA_TC_NEIGHBORHOOD, /*tsb.dimension()*/m_intrinsic_dim))
      : std::numeric_limits<std::size_t>::max();
#endif

    // Insert points until we find a point which is outside "star sphere"
    for (INS_iterator nn_it = ins_range.begin() ;
      nn_it != ins_range.end() ;
      ++nn_it)
    {
#ifdef CGAL_ALPHA_TC  
      ++number_of_attempts_to_insert_points;
      /*if (number_of_attempts_to_insert_points > MAX_NUM_INSERTED_POINTS)
      break;*/
#endif

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

        if (squared_star_sphere_radius_plus_margin
          && k_sqdist(center_pt, neighbor_pt)
             > *squared_star_sphere_radius_plus_margin)
             break;

        Tr_point proj_pt = project_point_and_compute_weight(
          neighbor_pt, neighbor_weight, tsb,
          local_tr_traits);

#ifdef CGAL_ALPHA_TC
        /*for (int i = 0 ; i < tsb.num_thickening_vectors() ; ++i)
        {
        FT c = coord(proj_pt, m_intrinsic_dim + i);
        if (c > 2.5*tsb.alpha_plus(i) + 2*m_half_sparsity  //CJTODO: il faut 2*tsb.alpha_plus(i)
        || c < 2.5*tsb.alpha_minus(i) - 2*m_half_sparsity) //CJTODO: il faut 2*tsb.alpha_minus(i)
        {
        goto end_of_insert_loop; // "continue" in n-2 for-loop
        }

        }*/
        /*if (max_sqdist_to_tangent_space != FT(0)
        && center_to_nbor_sqdist > max_sqdist_to_tangent_space)
        break;*/
#endif


#ifdef CGAL_TC_VERY_VERBOSE
        ++num_attempts_to_insert_points;
#endif


        Tr_vertex_handle vh = triangulation.insert_if_in_star(proj_pt, center_vertex);
        //Tr_vertex_handle vh = triangulation.insert(proj_pt);
        if (vh != Tr_vertex_handle())
        {
#ifdef CGAL_TC_VERY_VERBOSE
          ++num_inserted_points;
#endif
          if (verbose)
            std::cerr << "* Inserted point #" << neighbor_point_idx << "\n";

          vh->data() = neighbor_point_idx;

          // Let's recompute squared_star_sphere_radius_plus_margin
          if (triangulation.current_dimension() >= tangent_space_dim)
          {
            squared_star_sphere_radius_plus_margin = boost::none;
            // Get the incident cells and look for the biggest circumsphere
            std::vector<Tr_full_cell_handle> incident_cells;
            triangulation.incident_full_cells(
              center_vertex,
              std::back_inserter(incident_cells));
            for (typename std::vector<Tr_full_cell_handle>::iterator cit =
              incident_cells.begin(); cit != incident_cells.end(); ++cit)
            {
              Tr_full_cell_handle cell = *cit;
              if (triangulation.is_infinite(cell))
              {
                squared_star_sphere_radius_plus_margin = boost::none;
                break;
              }
              else
              {
                // Note that this uses the perturbed point since it uses
                // the points of the local triangulation
                Tr_point c = power_center(
                  boost::make_transform_iterator(
                  cell->vertices_begin(),
                  vertex_handle_to_point<Tr_point, Tr_vertex_handle>),
                  boost::make_transform_iterator(
                  cell->vertices_end(),
                  vertex_handle_to_point<Tr_point, Tr_vertex_handle>));

                FT sq_power_sphere_diam = 4 * point_weight(c);

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
              // "4*m_sq_half_sparsity" because both points can be perturbed
              squared_star_sphere_radius_plus_margin =
                *squared_star_sphere_radius_plus_margin + 4 * m_sq_half_sparsity;
#else
              // "2*m_half_sparsity" because both points can be perturbed
              squared_star_sphere_radius_plus_margin = CGAL::square(
                CGAL::sqrt(*squared_star_sphere_radius_plus_margin)
                + 2 * m_half_sparsity);
#endif
              // Save it in `m_squared_star_spheres_radii_incl_margin`
              m_squared_star_spheres_radii_incl_margin[i] = 
                *squared_star_sphere_radius_plus_margin;
            }
            else
            {
              m_squared_star_spheres_radii_incl_margin[i] = FT(-1);
            }
          }
        }
      }
    }

    return center_vertex;
  }

  void refresh_tangent_triangulation(
    std::size_t i, Points_ds const& updated_pts_ds, bool verbose = false)
  {
    if (verbose)
      std::cerr << "** Refreshing tangent tri #" << i << " **\n";

    if (m_squared_star_spheres_radii_incl_margin[i] == FT(-1))
      return compute_tangent_triangulation(i, verbose);
    
    Point center_point = compute_perturbed_point(i);
    // Among updated point, what is the closer from our center point?
    std::size_t closest_pt_index = 
      updated_pts_ds.query_ANN(center_point, 1, false).begin()->first;

    typename K::Construct_weighted_point_d k_constr_wp =
      m_k.construct_weighted_point_d_object();
    typename K::Power_distance_d k_power_dist = m_k.power_distance_d_object();
    
    // Construct a weighted point equivalent to the star sphere
    Weighted_point star_sphere = k_constr_wp(
      compute_perturbed_point(i),
      m_squared_star_spheres_radii_incl_margin[i]);
    Weighted_point closest_updated_point =
      compute_perturbed_weighted_point(closest_pt_index);

    // Is the "closest point" inside our star sphere?
    if (k_power_dist(star_sphere, closest_updated_point) <= FT(0))
      compute_tangent_triangulation(i, verbose);
  }

  void compute_tangent_triangulation(std::size_t i, bool verbose = false)
  {
    if (verbose)
      std::cerr << "** Computing tangent tri #" << i << " **\n";
    //std::cerr << "***********************************************\n";

    // No need to lock the mutex here since this will not be called while
    // other threads are perturbing the positions
    const Point center_pt = compute_perturbed_point(i);
    Tangent_space_basis &tsb = m_tangent_spaces[i];
    
#if defined(CGAL_TC_VERY_VERBOSE) && defined(CGAL_ALPHA_TC)
    std::cerr << "Base dimension, incl. thickening vectors: " 
      << tsb.dimension() << "\n";
#endif
    // Estimate the tangent space
    if (!m_are_tangent_spaces_computed[i])
    {
#ifdef CGAL_TC_EXPORT_NORMALS
      tsb = compute_tangent_space(center_pt, i, true/*normalize*/, &m_orth_spaces[i]);
#else
      tsb = compute_tangent_space(center_pt, i);
#endif
    }
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
    else if (m_perturb_tangent_space[i])
    {
#ifdef CGAL_TC_EXPORT_NORMALS
      tsb = compute_tangent_space(center_pt, i,
                                  true /*normalize_basis*/,
                                  &m_orth_spaces[i],
                                  true /*perturb*/);
#else
      tsb = compute_tangent_space(center_pt, i,
                                  true /*normalize_basis*/,
                                  NULL /*ortho basis*/,
                                  true /*perturb*/);
#endif
      m_perturb_tangent_space[i] = false;
    }
#endif

#if defined(CGAL_TC_PROFILING) && defined(CGAL_TC_VERY_VERBOSE)
    Wall_clock_timer t;
#endif
    int tangent_space_dim = tangent_basis_dim(i);
    Triangulation &local_tr =
      m_triangulations[i].construct_triangulation(tangent_space_dim);

    m_triangulations[i].center_vertex() =
      compute_star(i, center_pt, tsb, local_tr, verbose);

#if defined(CGAL_TC_PROFILING) && defined(CGAL_TC_VERY_VERBOSE)
    std::cerr << "  - triangulation construction: " << t.elapsed() << " s.\n";
    t.reset();
#endif

#ifdef CGAL_TC_VERY_VERBOSE
    std::cerr << "Inserted " << num_inserted_points << " points / " 
      << num_attempts_to_insert_points << " attemps to compute the star\n";
#endif
#ifdef CGAL_ALPHA_TC
    if (tsb.num_thickening_vectors() == 0)
      update_star__no_thickening_vectors(i);
    else
    {
      update_star__with_thickening_vector(i);
    }
#else
    update_star__no_thickening_vectors(i);
#endif

#if defined(CGAL_TC_PROFILING) && defined(CGAL_TC_VERY_VERBOSE)
    std::cerr << "  - update_star: " << t.elapsed() << " s.\n";
#endif
  }

  // Updates m_stars[i] directly from m_triangulations[i]
  void update_star__no_thickening_vectors(std::size_t i)
  {
    //***************************************************
    // Update the associated star (in m_stars)
    //***************************************************
    Star &star = m_stars[i];
    star.clear();
    Triangulation &local_tr = m_triangulations[i].tr();
    Tr_vertex_handle center_vertex = m_triangulations[i].center_vertex();
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
  void update_star__with_thickening_vector(std::size_t i)
  {
    //***************************************************
    // Parse the faces of the star and add the ones that are in the
    // restriction to alpha-Tp
    // Update the associated star (in m_stars)
    //***************************************************
    
    Triangulation &local_tr = m_triangulations[i].tr();
    int triangulation_dim = local_tr.current_dimension();
    Tr_traits const& local_tr_traits = local_tr.geom_traits();
    Tr_vertex_handle center_vertex = m_triangulations[i].center_vertex();
    Tangent_space_basis const& tsb = m_tangent_spaces[i];

#ifdef CGAL_TC_PERFORM_EXTRA_CHECKS
    if (triangulation_dim != tangent_basis_dim(i))
      std::cerr << "WARNING in update_star__with_thickening_vector: the "
                   "dimension of the local triangulation is different from "
                   "the dimension of the tangent space.\n";
#endif

    Star &star = m_stars[i];
    star.clear();
    int cur_dim_plus_1 = triangulation_dim + 1;

    std::vector<Tr_full_cell_handle> incident_cells;
    local_tr.incident_full_cells(
      center_vertex, std::back_inserter(incident_cells));

    typedef std::set<Tr_vertex_handle> DT_face; // DT face without center vertex (i)
    typedef std::set<Tr_vertex_handle> Neighbor_vertices;
    typedef std::map<DT_face, Neighbor_vertices> DT_faces_and_neighbors;

    // Maps that associate a m-face F and the points of its m+1-cofaces
    // (except the points of F). Those points are called its "neighbors".
    // N.B.: each m-face contains 'i', so 'i' is not stored in the faces
    // N.B.2: faces_and_neighbors[0] => dim 1, faces_and_neighbors[1] => dim 2
    std::vector<DT_faces_and_neighbors> faces_and_neighbors;
    faces_and_neighbors.resize(triangulation_dim);

    // Fill faces_and_neighbors
    // Let's first take care of the maximal simplices (dim = triangulation_dim)
    typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
    typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end = incident_cells.end();
    // For each cell
    for ( ; it_c != it_c_end ; ++it_c)
    {
      DT_face face;
      for (int j = 0 ; j < cur_dim_plus_1 ; ++j)
      {
        Tr_vertex_handle vh = (*it_c)->vertex(j);
        // Skip infinite simplices
        if (vh == local_tr.infinite_vertex())
          goto next_face;
        if (vh->data() != i)
          face.insert(vh);
      }
      // No co-faces => no neighbors
      faces_and_neighbors[triangulation_dim-1][face] = Neighbor_vertices();
next_face:
      ;
    }
    // Then the D-m-faces...
    int current_dim = triangulation_dim - 1;
    while (current_dim > 0)
    {
      // Let's fill faces_and_neighbors[current_dim-1]
      // (stores the current_dim-faces)
      DT_faces_and_neighbors& cur_faces_and_nghb =
        faces_and_neighbors[current_dim-1];

      typedef DT_faces_and_neighbors::const_iterator FaN_it;
      // Parse m+1-faces
      for (FaN_it it_k_p1_face = faces_and_neighbors[current_dim].begin(),
                  it_k_p1_face_end = faces_and_neighbors[current_dim].end() ;
           it_k_p1_face != it_k_p1_face_end ; ++it_k_p1_face)
      {
        DT_face const& k_p1_face = it_k_p1_face->first;

        // Add each m-face to cur_faces_and_nghb
        std::size_t n = current_dim + 1; // Not +2 since 'i' is not stored
        std::vector<bool> booleans(n, false);
        std::fill(booleans.begin() + 1, booleans.end(), true);
        do
        {
          DT_face k_face;
          Tr_vertex_handle remaining_vertex;
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
    // Test if V intersects the thickened tangent space
    current_dim = triangulation_dim;
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
        std::vector<Tr_point> P;
        P.reserve(current_DT_face.size() + 1);
        for (DT_face::const_iterator it = current_DT_face.begin(),
          it_end = current_DT_face.end() ; it != it_end ; ++it)
        {
          P.push_back((*it)->point());
        }
        P.push_back(center_vertex->point());
        
        // Q: vertices which are common neighbors of all vertices of P
        std::vector<Tr_point> Q;
        P.reserve(curr_neighbors.size());
        for (Neighbor_vertices::const_iterator it = curr_neighbors.begin(),
          it_end = curr_neighbors.end() ; it != it_end ; ++it)
        {
          Q.push_back((*it)->point());
        }

        bool does_intersect =
          does_voronoi_face_and_tangent_subspace_intersect(
            triangulation_dim,
            P, 
            Q, 
            tsb, 
            local_tr_traits);
        if (does_intersect)
        {
          // Get the indices of the face's points
          Incident_simplex face;
          DT_face::const_iterator it_vh = current_DT_face.begin();
          DT_face::const_iterator it_vh_end = current_DT_face.end();
          for ( ; it_vh != it_vh_end ; ++it_vh)
            face.insert((*it_vh)->data());

          star.push_back(face);

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
  }
#endif //CGAL_ALPHA_TC


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
#ifdef CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_SPHERE_2

    double tt[2] = {p[1], -p[0]};
    Vector t(2, &tt[0], &tt[2]);

    // Normalize t1 and t2
    typename K::Squared_length_d sqlen      = m_k.squared_length_d_object();
    typename K::Scaled_vector_d  scaled_vec = m_k.scaled_vector_d_object();

    Tangent_space_basis ts(i);
    ts.reserve(m_intrinsic_dim);
    ts.push_back(scaled_vec(t, FT(1)/CGAL::sqrt(sqlen(t))));
    m_are_tangent_spaces_computed[i] = true;

    return ts;

#elif defined(CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_SPHERE_3)

    double tt1[3] = {-p[1] - p[2], p[0], p[0]};
    double tt2[3] = {p[1] * tt1[2] - p[2] * tt1[1],
                     p[2] * tt1[0] - p[0] * tt1[2],
                     p[0] * tt1[1] - p[1] * tt1[0]};
    Vector t1(3, &tt1[0], &tt1[3]);
    Vector t2(3, &tt2[0], &tt2[3]);

    // Normalize t1 and t2
    typename K::Squared_length_d sqlen      = m_k.squared_length_d_object();
    typename K::Scaled_vector_d  scaled_vec = m_k.scaled_vector_d_object();

    Tangent_space_basis ts(i);
    ts.reserve(m_intrinsic_dim);
    ts.push_back(scaled_vec(t1, FT(1)/CGAL::sqrt(sqlen(t1))));
    ts.push_back(scaled_vec(t2, FT(1)/CGAL::sqrt(sqlen(t2))));

    m_are_tangent_spaces_computed[i] = true;

    return ts;

#elif defined(CGAL_TC_COMPUTE_TANGENT_PLANES_FOR_TORUS_D)

    Tangent_space_basis ts(i);
    ts.reserve(m_intrinsic_dim);
    for (int dim = 0 ; dim < m_intrinsic_dim ; ++dim)
    {
      std::vector<FT> tt(m_ambient_dim, 0.);
      tt[2*dim] = -p[2*dim + 1];
      tt[2*dim + 1] = p[2*dim];
      Vector t(2*m_intrinsic_dim, tt.begin(), tt.end());
      ts.push_back(t);
    }

    m_are_tangent_spaces_computed[i] = true;
    
    //return compute_gram_schmidt_basis(ts, m_k);
    return ts;
    //******************************* PCA *************************************
    
#else

    unsigned int num_points_for_pca = static_cast<unsigned int>(
      std::pow(BASE_VALUE_FOR_PCA, m_intrinsic_dim));

    // Kernel functors
    typename K::Construct_vector_d      constr_vec =
      m_k.construct_vector_d_object();
    typename K::Compute_coordinate_d    coord =
      m_k.compute_coordinate_d_object();
    typename K::Squared_length_d        sqlen =
      m_k.squared_length_d_object();
    typename K::Scaled_vector_d         scaled_vec =
      m_k.scaled_vector_d_object();
    typename K::Scalar_product_d        scalar_pdct =
      m_k.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec =
      m_k.difference_of_vectors_d_object();
    //typename K::Translated_point_d      transl =
    //  m_k.translated_point_d_object();

#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
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
        //  points_for_pca[nn_it->first], m_translations[nn_it->first]);
        mat_points(j, i) = CGAL::to_double(coord(points_for_pca[nn_it->first], i));
#ifdef CGAL_TC_ADD_NOISE_TO_TANGENT_SPACE
        mat_points(j, i) += m_random_generator.get_double(
            -0.5*m_half_sparsity, 0.5*m_half_sparsity);
#endif
#ifdef CGAL_TC_PERTURB_TANGENT_SPACE
        if (perturb)
          mat_points(j, i) += m_random_generator.get_double(
            -0.5*m_half_sparsity, 0.5*m_half_sparsity);
#endif
      }
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    Tangent_space_basis tsb(i); // p = compute_perturbed_point(i) here

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    for (int j = m_ambient_dim - 1 ;
         j >= m_ambient_dim - m_intrinsic_dim ;
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
      p_orth_space_basis->set_origin(i);
      for (int j = m_ambient_dim - m_intrinsic_dim - 1 ;
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
#if defined(CGAL_ALPHA_TC) && defined(CGAL_USE_A_FIXED_ALPHA)
    // Add the orthogonal vectors as "thickening vectors"
    for (int j = m_ambient_dim - m_intrinsic_dim - 1 ;
          j >= 0 ;
          --j)
    {
      Vector v = constr_vec(m_ambient_dim,
                            eig.eigenvectors().col(j).data(),
                            eig.eigenvectors().col(j).data() + m_ambient_dim);
      tsb.add_thickening_vector(
        normalize_vector(v, m_k), -CGAL_TC_ALPHA_VALUE, CGAL_TC_ALPHA_VALUE);
    }
#endif

    m_are_tangent_spaces_computed[i] = true;

    //*************************************************************************

    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, FT(1)/sqrt(sqlen(n)));
    //std::cerr << "IP = " << scalar_pdct(n, ts[0]) << " & " << scalar_pdct(n, ts[1]) << "\n";

    return tsb;
    
#endif

    /*
    // Alternative code (to be used later)
    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, FT(1)/sqrt(sqlen(n)));
    //Vector t1(12., 15., 65.);
    //Vector t2(32., 5., 85.);
    //Tangent_space_basis ts;
    //ts.reserve(m_intrinsic_dim);
    //ts.push_back(diff_vec(t1, scaled_vec(n, scalar_pdct(t1, n))));
    //ts.push_back(diff_vec(t2, scaled_vec(n, scalar_pdct(t2, n))));
    //ts = compute_gram_schmidt_basis(ts, m_k);
    //return ts;
    */
  }

  // Compute the space tangent to a simplex (p1, p2, ... pn)
  // CJTODO: Improve this?
  // Basically, it takes all the neighbor points to p1, p2... pn and runs PCA
  // on it. Note that most points are duplicated.
  Tangent_space_basis compute_tangent_space(
    const Indexed_simplex &s, bool normalize_basis = true)
  {
    unsigned int num_points_for_pca = static_cast<unsigned int>(
      std::pow(BASE_VALUE_FOR_PCA, m_intrinsic_dim));

    // Kernel functors
    typename K::Construct_vector_d      constr_vec =
      m_k.construct_vector_d_object();
    typename K::Compute_coordinate_d    coord =
      m_k.compute_coordinate_d_object();
    typename K::Squared_length_d        sqlen =
      m_k.squared_length_d_object();
    typename K::Scaled_vector_d         scaled_vec =
      m_k.scaled_vector_d_object();
    typename K::Scalar_product_d        scalar_pdct =
      m_k.scalar_product_d_object();
    typename K::Difference_of_vectors_d diff_vec =
      m_k.difference_of_vectors_d_object();
    //typename K::Translated_point_d      transl =
    //  m_k.translated_point_d_object();

    // One row = one point
    Eigen::MatrixXd mat_points(s.size()*num_points_for_pca, m_ambient_dim);
    unsigned int current_row = 0;

    for (Indexed_simplex::const_iterator it_index = s.begin();
          it_index != s.end() ; ++it_index)
    {
      const Point &p = m_points[*it_index];
      
#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
      KNS_range kns_range = m_points_ds_for_tse.query_ANN(
        p, num_points_for_pca, false);
      const Points &points_for_pca = m_points_for_tse;
#else
      KNS_range kns_range = m_points_ds.query_ANN(p, num_points_for_pca, false);
      const Points &points_for_pca = m_points;
#endif

      KNS_iterator nn_it = kns_range.begin();
      for ( ;
           current_row < num_points_for_pca && nn_it != kns_range.end() ;
           ++current_row, ++nn_it)
      {
        for (int i = 0 ; i < m_ambient_dim ; ++i)
        {
          //const Point p = transl(
          //  points_for_pca[nn_it->first], m_translations[nn_it->first]);
          mat_points(current_row, i) = 
            CGAL::to_double(coord(points_for_pca[nn_it->first], i));
        }
      }
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    Tangent_space_basis tsb;

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    for (int j = m_ambient_dim - 1 ;
          j >= m_ambient_dim - m_intrinsic_dim ;
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

    return tsb;
  }

  // Returns the dimension of the ith local triangulation
  // This is particularly useful for the alpha-TC
  int tangent_basis_dim(std::size_t i) const
  {
    return m_tangent_spaces[i].dimension();
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
    typename K::Construct_weighted_point_d k_constr_wp =
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
    typename K::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename K::Scaled_vector_d k_scaled_vec =
      m_k.scaled_vector_d_object();
    typename Tr_traits::Compute_coordinate_d coord =
      tr_traits.compute_coordinate_d_object();

    Point global_point = compute_perturbed_point(tsb.origin());
    for (int i = 0 ; i < m_intrinsic_dim ; ++i)
      global_point = k_transl(global_point,
                              k_scaled_vec(tsb[i], coord(p, i)));

#ifdef CGAL_ALPHA_TC
    Tangent_space_basis::Thickening_vectors const& tv = tsb.thickening_vectors();
    for (int i = 0 ; i < tv.size() ; ++i)
    {
      global_point = k_transl(
        global_point,
        k_scaled_vec(tv[i].vec, coord(p, m_intrinsic_dim + i)));
    }
#endif
    return global_point;
  }

  // Project the point in the tangent space
  // Resulting point coords are expressed in tsb's space
  Tr_bare_point project_point(const Point &p,
                              const Tangent_space_basis &tsb) const
  {
    typename K::Scalar_product_d scalar_pdct =
      m_k.scalar_product_d_object();
    typename K::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();
    
    Vector v = diff_points(p, compute_perturbed_point(tsb.origin()));

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    coords.reserve(tsb.dimension());
    for (std::size_t i = 0 ; i < m_intrinsic_dim ; ++i)
    {
      // Local coords are given by the scalar product with the vectors of tsb
      FT coord = scalar_pdct(v, tsb[i]);
      coords.push_back(coord);
    }
        
#ifdef CGAL_ALPHA_TC
    Tangent_space_basis::Thickening_vectors const& tv = tsb.thickening_vectors();
    for (int i = 0 ; i < tv.size() ; ++i)
    {
      FT coord = scalar_pdct(v, tv[i].vec);
      coords.push_back(coord);
    }
#endif
    return Tr_bare_point(static_cast<int>(
      coords.size()), coords.begin(), coords.end());
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  // Resulting point coords are expressed in tsb's space
  Tr_point project_point_and_compute_weight(const Weighted_point &wp,
                                            const Tangent_space_basis &tsb,
                                            const Tr_traits &tr_traits) const
  {
    typename K::Point_drop_weight_d k_drop_w =
      m_k.point_drop_weight_d_object();
    typename K::Point_weight_d k_point_weight =
      m_k.point_weight_d_object();
    return project_point_and_compute_weight(
      k_drop_w(wp), k_point_weight(wp), tsb, tr_traits);
  }

  // Same as above, with slightly different parameters
  Tr_point project_point_and_compute_weight(const Point &p, const FT w,
                                            const Tangent_space_basis &tsb,
                                            const Tr_traits &tr_traits) const
  {
    const int point_dim = m_k.point_dimension_d_object()(p);

    typename K::Construct_point_d constr_pt =
      m_k.construct_point_d_object();
    typename K::Scalar_product_d scalar_pdct =
      m_k.scalar_product_d_object();
    typename K::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();
    typename K::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();
    typename K::Construct_cartesian_const_iterator_d ccci =
      m_k.construct_cartesian_const_iterator_d_object();

    Point origin = compute_perturbed_point(tsb.origin());
    Vector v = diff_points(p, origin);

    // Same dimension? Then the weight is 0
    bool same_dim = (point_dim == tsb.dimension());

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(ccci(origin), ccci(origin, 0));
    coords.reserve(tsb.dimension());
    for (std::size_t i = 0 ; i < tsb.dimension() ; ++i)
    {
      // Local coords are given by the scalar product with the vectors of tsb
      FT c = scalar_pdct(v, tsb[i]);
      coords.push_back(c);

      // p_proj += c * tsb[i]
      if (!same_dim)
        for (int j = 0 ; j < point_dim ; ++j)
          p_proj[j] += c * coord(tsb[i], j);
    }
    
#ifdef CGAL_ALPHA_TC
    Tangent_space_basis::Thickening_vectors const& tv = tsb.thickening_vectors();
    for (int i = 0 ; i < tv.size() ; ++i)
    {
      FT c = scalar_pdct(v, tv[i].vec);
      coords.push_back(c);
      
      // p_proj += c * tv[i].vec
      if (!same_dim)
        for (int j = 0 ; j < point_dim ; ++j)
          p_proj[j] += c * coord(tv[i].vec, j);
    }
#endif

    // Same dimension? Then the weight is 0
    FT sq_dist_to_proj_pt = 0;
    if (!same_dim)
    {
      Point projected_pt = constr_pt(point_dim, p_proj.begin(), p_proj.end());
      sq_dist_to_proj_pt = m_k.squared_distance_d_object()(p, projected_pt);
    }
    
    return tr_traits.construct_weighted_point_d_object()
    (
      constr_pt(static_cast<int>(coords.size()), coords.begin(), coords.end()),
      w - sq_dist_to_proj_pt
    );
  }
  
  // Project all the points in the tangent space
  template <typename Indexed_point_range>
  std::vector<Tr_point> project_points_and_compute_weights(
    const Indexed_point_range &point_indices,
    const Tangent_space_basis &tsb,
    const Tr_traits &tr_traits) const
  {
    std::vector<Tr_point> ret;
    for (typename Indexed_point_range::const_iterator 
      it = point_indices.begin(), it_end = point_indices.end();
      it != it_end ; ++it)
    {
      ret.push_back(project_point_and_compute_weight(
        compute_perturbed_weighted_point(*it), tsb, tr_traits));
    }
    return ret;
  }

  // A simplex here is a local tri's full cell handle
  bool is_simplex_consistent(Tr_full_cell_handle fch, int cur_dim) const
  {
    Indexed_simplex c;
    for (int i = 0 ; i < cur_dim + 1 ; ++i)
    {
      std::size_t data = fch->vertex(i)->data();
      c.insert(data);
    }
    return is_simplex_consistent(c);
  }
  
  // A simplex here is a list of point indices
  template <typename IndexRange>
  double compute_simplex_fatness(IndexRange const& simplex) const
  {
    // Kernel functors
    typename K::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();
    typename K::Squared_distance_d sqdist =
      m_k.squared_distance_d_object();
    typename K::Difference_of_points_d diff_pts =
      m_k.difference_of_points_d_object();
    
    typename Tr_traits::Difference_of_points_d tr_diff_pts =
      m_triangulations[0].tr().geom_traits().difference_of_points_d_object();

    std::vector<std::size_t> s(simplex.begin(), simplex.end());
    std::size_t simplex_dim = s.size() - 1;

    // Compute basis
    Tangent_space_basis basis(s[0]);
    for (int j = 0 ; j < simplex_dim  ; ++j)
    {
      Vector e = diff_pts(
        compute_perturbed_point(s[j+1]), compute_perturbed_point(s[0]));
      basis.push_back(e);
    }
    basis = compute_gram_schmidt_basis(basis, m_k);

    // Compute the volume of the simplex: determinant
    Eigen::MatrixXd m(simplex_dim, simplex_dim);
    for (int j = 0 ; j < simplex_dim ; ++j)
    {
      Tr_vector v_j = tr_diff_pts(
        project_point(compute_perturbed_point(s[j+1]), basis), 
        project_point(compute_perturbed_point(s[0]), basis));
      for (int i = 0 ; i < simplex_dim ; ++i)
        m(j, i) = CGAL::to_double(coord(v_j, i));
    }
    double volume = 
      std::abs(m.determinant()) 
      / boost::math::factorial<double>(simplex_dim);

    // Compute the longest edge of the simplex
    CGAL::Combination_enumerator<int> combi(2, 0, simplex_dim+1);
    FT max_sq_length = FT(0);
    for ( ; !combi.finished() ; ++combi)
    {
      FT sq_length = sqdist(
        compute_perturbed_point(s[combi[0]]), 
        compute_perturbed_point(s[combi[1]]));
      if (sq_length > max_sq_length)
        max_sq_length = sq_length;
    }

    return volume / std::pow(CGAL::sqrt(max_sq_length), simplex_dim);
  }

  // A simplex here is a list of point indices
  // CJTODO: improve it like the other "is_simplex_consistent" below
  bool is_simplex_consistent(Indexed_simplex const& simplex) const
  {
#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
    if (is_one_of_the_coord_far_from_origin(simplex, CGAL_TC_ALVAREZ_SURFACE_WINDOW, 2))
      return true;
#endif

    // Check if the simplex is in the stars of all its vertices
    Indexed_simplex::const_iterator it_point_idx = simplex.begin();
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
      Incident_simplex is_to_find = simplex;
      is_to_find.erase(point_idx);

      // For each cell
      if (std::find(star.begin(), star.end(), is_to_find) == star.end())
        return false;
    }

    return true;
  }
  
  // A simplex here is a list of point indices
  // "s" contains all the points of the simplex except "center_point"
  // This function returns the points whose star doesn't contain the simplex
  // N.B.: the function assumes that the simplex is contained in 
  //       star(center_point)
  template <typename OutputIterator> // value_type = std::size_t
  bool is_simplex_consistent(
    std::size_t center_point,
    Incident_simplex const& s, // without "center_point"
    OutputIterator points_whose_star_does_not_contain_s,
    bool check_also_in_non_maximal_faces = false) const
  {
    Indexed_simplex full_simplex = s;
    full_simplex.insert(center_point);

#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
    if (is_one_of_the_coord_far_from_origin(full_simplex, CGAL_TC_ALVAREZ_SURFACE_WINDOW, 2))
      return true;
#endif

    // Check if the simplex is in the stars of all its vertices
    Incident_simplex::const_iterator it_point_idx = s.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for ( ; it_point_idx != s.end() ; ++it_point_idx)
    {
      std::size_t point_idx = *it_point_idx;
      // Don't check infinite simplices
      if (point_idx == std::numeric_limits<std::size_t>::max())
        continue;

      Star const& star = m_stars[point_idx];

      // What we're looking for is full_simplex \ point_idx
      Incident_simplex is_to_find = full_simplex;
      is_to_find.erase(point_idx);

      if (check_also_in_non_maximal_faces)
      {
        // For each simplex "is" of the star, check if ic_to_simplex is
        // included in "is"
        bool found = false;
        for (Star::const_iterator is = star.begin(), is_end = star.end() ;
          !found && is != is_end ; ++is)
        {
          if (std::includes(is->begin(), is->end(), 
                            is_to_find.begin(), is_to_find.end()))
            found = true;
        }

        if (!found)
          *points_whose_star_does_not_contain_s++ = point_idx;
      }
      else
      {
        // Does the star contain is_to_find?
        if (std::find(star.begin(), star.end(), is_to_find) == star.end())
          *points_whose_star_does_not_contain_s++ = point_idx;
      }
    }

    return true;
  }

  // A simplex here is a list of point indices
  // It looks for s in star(p).
  // "s" contains all the points of the simplex except p.
  bool is_simplex_in_star(
    std::size_t p,
    Incident_simplex const& s,
    bool check_also_in_non_maximal_faces = true) const
  {
    Star const& star = m_stars[p];

    if (check_also_in_non_maximal_faces)
    {
      // For each simplex "is" of the star, check if ic_to_simplex is
      // included in "is"
      bool found = false;
      for (Star::const_iterator is = star.begin(), is_end = star.end() ;
        !found && is != is_end ; ++is)
      {
        if (std::includes(is->begin(), is->end(), s.begin(), s.end()))
          found = true;
      }

      return found;
    }
    else
    {
      return !(std::find(star.begin(), star.end(), s) == star.end());
    }
  }

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for try_to_solve_inconsistencies_in_a_local_triangulation function
  class Try_to_solve_inconsistencies_in_a_local_triangulation
  {
    Tangential_complex & m_tc;
    tbb::combinable<std::size_t> &m_num_inconsistencies;
    tbb::combinable<std::vector<std::size_t> > &m_updated_points;

  public:
    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
      Tangential_complex &tc,
      tbb::combinable<std::size_t> &num_inconsistencies,
      tbb::combinable<std::vector<std::size_t> > &updated_points)
    : m_tc(tc),
      m_num_inconsistencies(num_inconsistencies),
      m_updated_points(updated_points)
    {}

    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
      const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc),
      m_num_inconsistencies(ctt.m_num_inconsistencies),
      m_updated_points(ctt.m_updated_points)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
      {
        m_num_inconsistencies.local() +=
          m_tc.try_to_solve_inconsistencies_in_a_local_triangulation(
            i, std::back_inserter(m_updated_points.local()));
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
    typename K::Point_to_vector_d k_pt_to_vec =
      m_k.point_to_vector_d_object();
    CGAL::Random_points_in_ball_d<Point>
      tr_point_in_ball_generator(
        m_ambient_dim, m_half_sparsity);
    // Parallel
#  if defined(CGAL_LINKED_WITH_TBB) && defined(CGAL_TC_GLOBAL_REFRESH)
    Vector transl = k_pt_to_vec(*tr_point_in_ball_generator++);
    m_p_perturb_mutexes[point_idx].lock();
    m_translations[point_idx] = transl;
    m_p_perturb_mutexes[point_idx].unlock();
    // Sequential
#  else
    m_translations[point_idx] = k_pt_to_vec(*tr_point_in_ball_generator++);
#  endif

# else // CGAL_TC_PERTURB_POSITION_TANGENTIAL
    const Tr_traits &local_tr_traits =
      m_triangulations[point_idx].tr().geom_traits();
    typename Tr_traits::Compute_coordinate_d coord =
      local_tr_traits.compute_coordinate_d_object();
    typename K::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename K::Construct_vector_d k_constr_vec =
      m_k.construct_vector_d_object();
    typename K::Scaled_vector_d k_scaled_vec =
      m_k.scaled_vector_d_object();

    CGAL::Random_points_in_ball_d<Tr_bare_point>
      tr_point_in_ball_generator(
        m_intrinsic_dim, 
        m_random_generator.get_double(0., m_half_sparsity));

    Tr_point local_random_transl =
      local_tr_traits.construct_weighted_point_d_object()(
        *tr_point_in_ball_generator++, 0);
    Translation_for_perturb global_transl = k_constr_vec(m_ambient_dim);
    const Tangent_space_basis &tsb = m_tangent_spaces[point_idx];
    for (int i = 0 ; i < m_intrinsic_dim ; ++i)
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

  // Return true if inconsistencies were found
  template <typename OutputIt>
  bool try_to_solve_inconsistencies_in_a_local_triangulation(
    std::size_t tr_index, 
    OutputIt perturbed_pts_indices = CGAL::Emptyset_iterator())
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
      if (is_infinite(incident_simplex))
        continue;

      Indexed_simplex c = incident_simplex;
      c.insert(tr_index); // Add the missing index

//*****************************************************************************
// STRATEGY 1: perturb all the points of the first inconsistent simplex
//*****************************************************************************
#ifdef CGAL_TC_PERTURB_THE_SIMPLEX_ONLY
      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        for (Indexed_simplex::const_iterator it = c.begin();
             it != c.end() ; ++it)
        {
          perturb(*it);
          *perturbed_pts_indices++ = *it;
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
        *perturbed_pts_indices++ = idx;

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
          *perturbed_pts_indices++ = *it;
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
          CGAL_TC_NUMBER_OF_PERTURBED_POINTS(m_intrinsic_dim));
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
          *perturbed_pts_indices++ = *it;
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
        {
          perturb(tr_index);
          *perturbed_pts_indices++ = tr_index;
        }
        else
        {
          Indexed_simplex::const_iterator it_idx = c.begin();
          std::advance(it_idx, rnd - 1);
          perturb(*it_idx);
          *perturbed_pts_indices++ = *it_idx;
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


  bool filter_inconsistent_simplices_in_a_local_triangulation(
                                                          std::size_t tr_index)
  {
    bool is_inconsistent = false;

    Star &star = m_stars[tr_index];
    Triangulation const& tr    = m_triangulations[tr_index].tr();
    Tr_vertex_handle center_vh = m_triangulations[tr_index].center_vertex();
    const Tr_traits &local_tr_traits = tr.geom_traits();
    int cur_dim = tr.current_dimension();

    // For each incident simplex
    Star::iterator it_inc_simplex = star.begin();
    while (it_inc_simplex != star.end())
    {
      const Incident_simplex &incident_simplex = *it_inc_simplex;

      // Don't check infinite cells
      if (is_infinite(incident_simplex))
      {
        ++it_inc_simplex;
        continue;
      }

      Indexed_simplex c = incident_simplex;
      c.insert(tr_index); // Add the missing index

      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        // Compute T, the space tangent to "c"
        Tangent_space_basis tsb = compute_tangent_space(c);

        std::size_t center_index = *(c.begin());
        // Compute the star of "center_index" in the RT restricted to T
        tsb.set_origin(center_index);
        Triangulation tr(m_intrinsic_dim);
        Tr_vertex_handle center_vh = compute_star(
          center_index, compute_perturbed_point(center_index), tsb, tr);

        // Get incident cells
        std::vector<Tr_full_cell_handle> incident_cells;
        tr.incident_full_cells(
          center_vh, std::back_inserter(incident_cells));

        // Check if "c" is among them
        typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
        typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end = incident_cells.end();
        bool found = false;
        // For each cell
        for (; !found && it_c != it_c_end ; ++it_c)
        {
          Incident_simplex incident_simplex;
          for (int j = 0 ; j < tr.current_dimension() + 1 ; ++j)
          {
            std::size_t index = (*it_c)->vertex(j)->data();
            incident_simplex.insert(index);
          }
          if (incident_simplex == c)
            found = true;
        }
        // If the simplex is not in the new star, we remove it
        if (!found)
        {
          if (it_inc_simplex != star.end() - 1)
          {
            *it_inc_simplex = star.back();
            star.pop_back();
          }
          else
          {
            star.pop_back();
            it_inc_simplex = star.end(); // end the loop
          }
        }
        else
          ++it_inc_simplex;
      }
      // Consistent
      else
        ++it_inc_simplex;
    }

    return is_inconsistent;
  }

  // 1st line: number of points
  // Then one point per line
  std::ostream &export_point_set(
    std::ostream & os,
    bool use_perturbed_points = false,
    const char *coord_separator = " ") const
  {
    if (use_perturbed_points)
    {
      std::vector<Point> perturbed_points;
      perturbed_points.reserve(m_points.size());
      for (std::size_t i = 0 ; i < m_points.size() ; ++i)
        perturbed_points.push_back(compute_perturbed_point(i));

      return export_point_set(
        m_k, perturbed_points, os, coord_separator);
    }
    else
    {
      return export_point_set(
        m_k, m_points, os, coord_separator);
    }
  }

  template<typename ProjectionFunctor = CGAL::Identity<Point> >
  std::ostream &export_vertices_to_off(
    std::ostream & os, std::size_t &num_vertices,
    bool use_perturbed_points = false, 
    ProjectionFunctor const& point_projection = ProjectionFunctor()) const
  {
    if (m_points.empty())
    {
      num_vertices = 0;
      return os;
    }

    // If m_intrinsic_dim = 1, we output each point two times
    // to be able to export each segment as a flat triangle with 3 different
    // indices (otherwise, Meshlab detects degenerated simplices)
    const int N = (m_intrinsic_dim == 1 ? 2 : 1);

    // Kernel functors
    typename K::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();

#ifdef CGAL_TC_EXPORT_ALL_COORDS_IN_OFF
    int num_coords = m_ambient_dim;
#else
    int num_coords = min(m_ambient_dim, 3);
#endif

#ifdef CGAL_TC_EXPORT_NORMALS
    OS_container::const_iterator it_os = m_orth_spaces.begin();
#endif
    typename Points::const_iterator it_p = m_points.begin();
    typename Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for (std::size_t i = 0 ; it_p != it_p_end ; ++it_p, ++i)
    {
      Point p = point_projection(
        use_perturbed_points ? compute_perturbed_point(i) : *it_p);
      for (int ii = 0 ; ii < N ; ++ii)
      {
        int j = 0;
#if CGAL_TC_BETTER_EXPORT_FOR_FLAT_TORUS
        // For flat torus
        os << (2 + 1 * CGAL::to_double(coord(p, 0))) * CGAL::to_double(coord(p, 2)) << " "
           << (2 + 1 * CGAL::to_double(coord(p, 0))) * CGAL::to_double(coord(p, 3)) << " "
           << 1 * CGAL::to_double(coord(p, 1));
#else
        for ( ; j < num_coords ; ++j)
          os << CGAL::to_double(coord(p, j)) << " ";
#endif
        if (j == 2)
          os << "0";

#ifdef CGAL_TC_EXPORT_NORMALS
        for (j = 0 ; j < num_coords ; ++j)
          os << " " << CGAL::to_double(coord(*it_os->begin(), j));
#endif
        os << "\n";
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
    const Indexed_simplex &simplex)
  {
    Incident_simplex incident_simplex = simplex;
    incident_simplex.erase(index); // Remove the center index

    Star &star = m_stars[index];

    Indexed_simplex::const_iterator it_point_idx = simplex.begin();
    Indexed_simplex::const_iterator it_point_idx_end = simplex.end();
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
    const Indexed_simplex &inconsistent_simplex)
  {
    CGAL_assertion_code(
      Indexed_simplex inc_s_minus_p = inconsistent_simplex;
      inc_s_minus_p.erase(p_idx);
      Indexed_simplex inc_s_minus_q = inconsistent_simplex;
      inc_s_minus_q.erase(q_idx);
    );
    CGAL_assertion(std::find(m_stars[p_idx].begin(), m_stars[p_idx].end(),
                             inc_s_minus_p) != m_stars[p_idx].end());
    CGAL_assertion(std::find(m_stars[q_idx].begin(), m_stars[q_idx].end(),
                             inc_s_minus_q) == m_stars[q_idx].end());

    typename K::Point_drop_weight_d k_drop_w =
      m_k.point_drop_weight_d_object();
    typename K::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename K::Squared_distance_d k_sqdist =
      m_k.squared_distance_d_object();
    typename K::Difference_of_points_d k_diff_pts =
      m_k.difference_of_points_d_object();
    typename K::Scalar_product_d k_scalar_pdct =
      m_k.scalar_product_d_object();
    typename K::Construct_weighted_point_d k_constr_wp =
      m_k.construct_weighted_point_d_object();
    typename K::Power_distance_d k_power_dist =
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

    Indexed_simplex::const_iterator it_point_idx =
                                                  inconsistent_simplex.begin();
    Indexed_simplex::const_iterator it_point_idx_end =
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

    // CJTODO DEBUG ====================
    /*{
    INS_range ins_range = m_points_ds.query_incremental_ANN(k_drop_w(global_Cp));
    for (INS_iterator nn_it = ins_range.begin() ;
         nn_it != ins_range.end() ;
         ++nn_it)
    {
      FT neighbor_sqdist = nn_it->second;

      // When we're sure we got all the potential points, break
      //if (neighbor_sqdist > circumsphere_sqradius_p + m_sq_half_sparsity)
      //  break;

      std::size_t neighbor_point_idx = nn_it->first;
      FT point_to_Cp_power_sqdist = k_power_dist(
        global_Cp, compute_perturbed_weighted_point(neighbor_point_idx));
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
      }*/
    }
    CGAL_assertion_msg(!inside_pt_indices.empty(),
      "There should be at least one vertex inside the sphere");

    // CJTODO DEBUG
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
        std::cerr << "\n";
      }
      std::cerr << "\n";
    }*/

    // CJTODO DEBUG
    // If co-intrinsic dimension = 1, let's compare normals
    /*if (m_ambient_dim - m_intrinsic_dim == 1)
    {
      typename K::Scaled_vector_d k_scaled_vec =
        m_k.scaled_vector_d_object();
      typename K::Squared_length_d k_sqlen =
        m_k.squared_length_d_object();
      Vector pq = k_diff_pts(
        compute_perturbed_point(q_idx), compute_perturbed_point(p_idx));
      pq = k_scaled_vec(pq, FT(1)/sqrt(k_sqlen(pq)));
      FT dot_product_1 = std::abs(
          k_scalar_pdct(m_orth_spaces[p_idx][0], pq));
      FT dot_product_2 = std::abs(
          k_scalar_pdct(m_orth_spaces[q_idx][0], pq));
      csv_stream << inside_pt_indices.size() << " ; ";
      csv_stream << dot_product_1 << " ; " << dot_product_2;
      csv_stream << "\n";
    }*/
    
    // CJTODO DEBUG
    if (inside_pt_indices.size() > 1)
    {
      std::cerr << "Warning: " << inside_pt_indices.size() << " insiders in "
        << inconsistent_simplex.size() - 1 << " simplex\n";
      
      // If co-intrinsic dimension = 1, let's compare normals
      /*if (m_ambient_dim - m_intrinsic_dim == 1)
      {
        std::cerr << "(dot product between normals = ";
        Indexed_simplex::const_iterator it_v = 
          inconsistent_simplex.begin();
        std::size_t i1 = *it_v;
        ++it_v;
        for ( ; it_v != inconsistent_simplex.end() ; ++it_v)
        {
          FT dot_products_between_normals =
            k_scalar_pdct(m_tangent_spaces[i1][0], m_tangent_spaces[*it_v][0]);
          std::cerr << dot_products_between_normals << ", ";
          //csv_stream << " ; " <<dot_products_between_normals;
        }
        std::cerr << "\n";
        //csv_stream << "\n";
      }*/
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
          (FT(2)*k_scalar_pdct(k_diff_pts(cq, cp), k_diff_pts(ti, pt_p)));
#else
        FT a =
          (k_sqdist(cp, ti) - k_sqdist(cp, pt_p)) /
          (FT(2)*k_scalar_pdct(k_diff_pts(cq, cp), k_diff_pts(ti, pt_p)));
#endif

        if (a < min_a)
        {
          min_a = a;
          inside_point_idx = idx;
        }
      }

      // CJTODO DEBUG ====================
      /*{
      typename K::Scaled_vector_d scaled_vec = m_k.scaled_vector_d_object();
      typename K::Point_weight_d k_weight = m_k.point_weight_d_object();
      Weighted_point C = k_constr_wp(
        k_transl(k_drop_w(global_Cp), scaled_vec(k_diff_pts(k_drop_w(global_Cq), k_drop_w(global_Cp)), min_a)),
        k_sqdist(k_transl(k_drop_w(global_Cp), scaled_vec(k_diff_pts(k_drop_w(global_Cq), k_drop_w(global_Cp)), min_a)), pt_p));
      INS_range ins_range = m_points_ds.query_incremental_ANN(k_drop_w(C));
      for (INS_iterator nn_it = ins_range.begin() ;
           nn_it != ins_range.end() ;
           ++nn_it)
      {
        FT neighbor_sqdist = nn_it->second;

        // When we're sure we got all the potential points, break
        if (neighbor_sqdist > k_weight(C) + m_sq_half_sparsity)
          break;

        std::size_t neighbor_point_idx = nn_it->first;
        FT point_to_C_power_sqdist =
          k_power_dist(C, compute_perturbed_weighted_point(neighbor_point_idx));
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
    Indexed_simplex new_simplex = inconsistent_simplex;
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
    std::size_t tr_index, const Incident_simplex &incident_simplex)
  {
    bool inconsistencies_found = false;

    // Don't check infinite simplices
    if (is_infinite(incident_simplex))
      return false;

    Indexed_simplex simplex = incident_simplex;
    simplex.insert(tr_index);

    // Check if the simplex is in the stars of all its vertices
    Incident_simplex::const_iterator it_point_idx = incident_simplex.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for ( ; it_point_idx != incident_simplex.end() ; ++it_point_idx)
    {
      std::size_t point_idx = *it_point_idx;

      Star const& star = m_stars[point_idx];

      // What we're looking for is "simplex" \ point_idx
      Incident_simplex is_to_find = simplex;
      is_to_find.erase(point_idx);

      if (std::find(star.begin(), star.end(), is_to_find) == star.end())
      {
        solve_inconsistency_by_adding_higher_dimensional_simplices(
          tr_index, *it_point_idx, simplex);
        inconsistencies_found = true;
        break;
      }
      // CJTODO DEBUG
      /*else if (m_ambient_dim - m_intrinsic_dim == 1)
      {
        typename K::Difference_of_points_d k_diff_pts =
          m_k.difference_of_points_d_object();
        typename K::Scaled_vector_d k_scaled_vec =
          m_k.scaled_vector_d_object();
        typename K::Squared_length_d k_sqlen =
          m_k.squared_length_d_object();
        typename K::Scalar_product_d k_scalar_pdct =
          m_k.scalar_product_d_object();
        Vector pq = k_diff_pts(
          compute_perturbed_point(*it_point_idx), compute_perturbed_point(tr_index));
        pq = k_scaled_vec(pq, FT(1)/sqrt(k_sqlen(pq)));
        FT dot_product_1 = std::abs(
            k_scalar_pdct(m_orth_spaces[tr_index][0], pq));
        FT dot_product_2 = std::abs(
            k_scalar_pdct(m_orth_spaces[*it_point_idx][0], pq));
        csv_stream << "0 ; ";
        csv_stream << dot_product_1 << " ; " << dot_product_2;
        csv_stream << "\n";
      }*/
    }

    return inconsistencies_found;
  }

  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  // Note the computation is made in local coordinates. "tsb"'s vectors are not
  // used because these vectors become (0..., 1..., 0) in local coordinates.
  template <typename Weighted_point_range_a, typename Weighted_point_range_b>
  CGAL::Quadratic_program_solution<ET> 
    compute_voronoi_face_and_tangent_subspace_LP_problem(
    int points_dim,
    Weighted_point_range_a const& P,
    Weighted_point_range_b const& Q,
    Tangent_space_basis const& tsb,
    const Tr_traits &tr_traits) const
  {
    // Notations:
    // Fv: Voronoi k-face
    // Fd: dual, (D-k)-face of Delaunay (p0, p1, ..., pn)
    
    typename Tr_traits::Point_drop_weight_d drop_w =
      tr_traits.point_drop_weight_d_object();
    typename Tr_traits::Point_weight_d point_weight =
      tr_traits.point_weight_d_object();
    typename Tr_traits::Scalar_product_d scalar_pdct = 
      tr_traits.scalar_product_d_object();
    typename Tr_traits::Point_to_vector_d pt_to_vec = 
      tr_traits.point_to_vector_d_object();
    typename Tr_traits::Compute_coordinate_d coord = 
      tr_traits.compute_coordinate_d_object();
    
    typename K::Compute_coordinate_d k_coord = 
      m_k.compute_coordinate_d_object();

    std::size_t card_P = P.size();
    std::size_t card_Q = Q.size();

    // Linear solver
    typedef CGAL::Quadratic_program<FT> Linear_program;
    typedef CGAL::Quadratic_program_solution<ET> LP_solution;

    Linear_program lp(CGAL::SMALLER, false);
    int current_row = 0;

    //=========== First set of equations ===========
    // For points pi in P
    //   2(p0 - pi).x = p0^2 - wght(p0) - pi^2 + wght(pi)
    typename Weighted_point_range_a::const_iterator it_p = P.begin();
    Tr_point const& p0 = *it_p;
    FT const w0 = point_weight(p0);
    FT p0_dot_p0 = scalar_pdct(pt_to_vec(drop_w(p0)), pt_to_vec(drop_w(p0)));
    ++it_p;
    for (typename Weighted_point_range_a::const_iterator it_p_end = P.end() ;
         it_p != it_p_end ; ++it_p)
    {
      Tr_point const& pi = *it_p;
      FT const wi = point_weight(pi);

      for (int k = 0 ; k < points_dim ; ++k)
        lp.set_a(k, current_row, 2*(coord(p0, k) - coord(pi, k)));

      FT pi_dot_pi = scalar_pdct(pt_to_vec(drop_w(pi)), pt_to_vec(drop_w(pi)));
      lp.set_b(current_row, p0_dot_p0 - pi_dot_pi - w0 + wi);
      lp.set_r(current_row, CGAL::EQUAL);

      ++current_row;
    }

    //=========== Second set of equations ===========
    // For each point qi in Q
    //  2(qi - p0).x <= qi^2 - wght(pi) - p0^2 + wght(p0)
    for (typename Weighted_point_range_b::const_iterator it_q = Q.begin(),
                                             it_q_end = Q.end() ;
         it_q != it_q_end ; ++it_q)
    {
      Tr_point const& qi = *it_q;
      FT const wi = point_weight(qi);

      for (int k = 0 ; k < points_dim ; ++k)
        lp.set_a(k, current_row, 2*(coord(qi, k) - coord(p0, k)));

      FT qi_dot_qi = scalar_pdct(pt_to_vec(drop_w(qi)), pt_to_vec(drop_w(qi)));
      lp.set_b(current_row, qi_dot_qi - wi - p0_dot_p0 + w0);

      ++current_row;
    }

    //=========== Third set of equations ===========
    // For each thickening vector bj of TSB, 
    // x.bj <= alpha_plus and >= alpha_minus
    // where bj is in the TSB => bj = (0..., 1..., 0) (1 is at the jth position)
    //   x.bj  <=  alpha_plus
    //   -x.bj <= -alpha_minus
    std::size_t j = points_dim - tsb.num_thickening_vectors();
    for (Tangent_space_basis::Thickening_vectors::const_iterator 
           it_tv = tsb.thickening_vectors().begin(),
           it_tv_end = tsb.thickening_vectors().end() ;
         it_tv != it_tv_end ; ++it_tv)
    {
      Tangent_space_basis::Thickening_vector const& bj = *it_tv;

      for (int k = 0 ; k < points_dim ; ++k)
      {
        lp.set_a(k, current_row    , (j == k ? 1. : 0.));
        lp.set_a(k, current_row + 1, (j == k ? -1. : 0.));
      }

      lp.set_b(current_row    ,  bj.alpha_plus);
      lp.set_b(current_row + 1, -bj.alpha_minus);

      current_row += 2;
      ++j;
    }

    //=========== Other LP parameters ===========
    lp.set_c(0, 1); // Minimize x[0]

    //=========== Solve =========================
    LP_solution solution = CGAL::solve_linear_program(lp, ET());
    return solution;
  }

  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  template <typename Weighted_point_range_a, typename Weighted_point_range_b>
  bool does_voronoi_face_and_tangent_subspace_intersect(
    int points_dim,
    Weighted_point_range_a const& P,
    Weighted_point_range_b const& Q,
    Tangent_space_basis const& tsb,
    const Tr_traits &tr_traits) const
  {
    return compute_voronoi_face_and_tangent_subspace_LP_problem(
      points_dim, P, Q, tsb, tr_traits).status() == CGAL::QP_OPTIMAL;
  }
  
  // Returns any point of the intersection between aff(voronoi_cell) and a
  // tangent space.
  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Return value: the point coordinates are expressed in the tsb base
  template <typename Weighted_point_range>
  boost::optional<Tr_bare_point> 
  compute_aff_of_voronoi_face_and_tangent_subspace_intersection(
    int points_dim,
    Weighted_point_range const& P,
    Tangent_space_basis const& tsb,
    const Tr_traits &tr_traits) const
  {
    // As we're only interested by aff(v), Q is empty
    return compute_voronoi_face_and_tangent_subspace_intersection(
      points_dim, P, std::vector<typename Weighted_point_range::value_type>(), 
      tsb, tr_traits);
  }
  
  // Returns any point of the intersection between a Voronoi cell and a
  // tangent space.
  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  // Return value: the point coordinates are expressed in the tsb base
  template <typename Weighted_point_range_a, typename Weighted_point_range_b>
  boost::optional<Tr_bare_point> 
  compute_voronoi_face_and_tangent_subspace_intersection(
    int points_dim,
    Weighted_point_range_a const& P,
    Weighted_point_range_b const& Q,
    Tangent_space_basis const& tsb,
    const Tr_traits &tr_traits) const
  {
    typedef CGAL::Quadratic_program_solution<ET> LP_solution;
    
    LP_solution sol = compute_voronoi_face_and_tangent_subspace_LP_problem(
      points_dim, P, Q, tsb, tr_traits);

    boost::optional<Tr_bare_point> ret;
    if (sol.status() == CGAL::QP_OPTIMAL)
    {
      std::vector<FT> p;
      p.reserve(points_dim);
      for (LP_solution::Variable_value_iterator 
        it_v = sol.variable_values_begin(),
        it_v_end = sol.variable_values_end() ;
        it_v != it_v_end ; ++it_v)
      {
        p.push_back(to_double(*it_v));
      }
      ret = tr_traits.construct_point_d_object()(points_dim, p.begin(), p.end());
    }
    else
    {
      ret = boost::none;
    }

    return ret;
  }

  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  template <typename Indexed_point_range_a, typename Indexed_point_range_b>
  bool does_voronoi_face_and_fixed_alpha_tangent_subspace_intersect(
      std::size_t center_pt_index,
      Indexed_point_range_a const& P,
      Indexed_point_range_b const& Q,
      Orthogonal_space_basis const& orthogonal_subspace_basis,
      FT alpha) const
  {
    // Notations:
    // Fv: Voronoi k-face
    // Fd: dual, (D-k)-face of Delaunay (p0, p1, ..., pn)

    typename K::Scalar_product_d scalar_pdct = m_k.scalar_product_d_object();
    typename K::Point_to_vector_d pt_to_vec = m_k.point_to_vector_d_object();
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

    Point center_pt = compute_perturbed_point(center_pt_index);
    int const ambient_dim = m_k.point_dimension_d_object()(center_pt);

    std::size_t card_P = P.size();
    std::size_t card_Q = Q.size();

    // Linear solver
    typedef CGAL::Quadratic_program<FT> Linear_program;
    typedef CGAL::Quadratic_program_solution<ET> LP_solution;

    Linear_program lp(CGAL::SMALLER, false);
    int current_row = 0;

    //=========== First set of equations ===========
    // For points pi in P
    //   2(p0 - pi).x = p0^2 - w0 - pi^2 + wi
    Point const& p0 = center_pt;
    FT const w0 = m_weights[center_pt_index];
    FT p0_dot_p0 = scalar_pdct(pt_to_vec(p0), pt_to_vec(p0));

    for (typename Indexed_point_range_a::const_iterator it_p = P.begin(),
                                                        it_p_end = P.end() ;
         it_p != it_p_end ; ++it_p)
    {
      Point pi;
      FT wi;
      compute_perturbed_weighted_point(*it_p, pi, wi);

      for (int k = 0 ; k < ambient_dim ; ++k)
        lp.set_a(k, current_row, 2*(coord(p0, k) - coord(pi, k)));

      FT pi_dot_pi = scalar_pdct(pt_to_vec(pi), pt_to_vec(pi));
      lp.set_b(current_row, p0_dot_p0 - pi_dot_pi - w0 + wi);
      lp.set_r(current_row, CGAL::EQUAL);

      ++current_row;
    }

    // CJTODO: this code might be useful for Option 1
    /*CGAL::Combination_enumerator<int> pi_pj(2, 0, static_cast<int>(card_P));
    for ( ; !pi_pj.finished() ; ++pi_pj)
    {
      Point const& pi = P[pi_pj[0]];
      FT wi = all_weights[pi_pj[0]];
      Point const& pj = P[pi_pj[1]];
      FT wj = all_weights[pi_pj[1]];

      for (int k = 0 ; k < ambient_dim ; ++k)
      {
        FT a = 2*(coord(pi, k) + coord(pj, k));
        lp.set_a(k, current_row    , -a);
        lp.set_a(k, current_row + 1,  a);
      }

      FT b = scalar_pdct(pi, pi) - wi - scalar_pdct(pj, pj) + wj;
      lp.set_b(current_row    , -b);
      lp.set_b(current_row + 1,  b);

      current_row += 2;
    }*/

    //=========== Second set of equations ===========
    // For each point qi in Q
    //  2(qi - p0).x <= qi^2 - wi - p0^2 + w0
    for (typename Indexed_point_range_b::const_iterator it_q = Q.begin(),
                                                        it_q_end = Q.end() ;
         it_q != it_q_end ; ++it_q)
    {
      Point qi;
      FT wi;
      compute_perturbed_weighted_point(*it_q, qi, wi);

      for (int k = 0 ; k < ambient_dim ; ++k)
        lp.set_a(k, current_row, 2*(coord(qi, k) - coord(p0, k)));

      FT qi_dot_qi = scalar_pdct(pt_to_vec(qi), pt_to_vec(qi));
      lp.set_b(current_row, qi_dot_qi - wi - p0_dot_p0 + w0);

      ++current_row;
    }

    //=========== Third set of equations ===========
    // For each vector bi of OSB, (x-p).bi <= alpha and >= -alpha
    // p is the origin of the basis
    //     bi.x <=  bi.p + alpha
    //    -bi.x <= -bi.p + alpha
    for (Orthogonal_space_basis::const_iterator it_osb =
          orthogonal_subspace_basis.begin(),
          it_osb_end = orthogonal_subspace_basis.end() ;
         it_osb != it_osb_end ; ++it_osb)
    {
      Vector const& bi = *it_osb;

      for (int k = 0 ; k < ambient_dim ; ++k)
      {
        lp.set_a(k, current_row    ,  coord(bi, k));
        lp.set_a(k, current_row + 1, -coord(bi, k));
      }

      FT bi_dot_p = scalar_pdct(bi, 
        pt_to_vec(compute_perturbed_point(orthogonal_subspace_basis.origin())));
      lp.set_b(current_row    ,  bi_dot_p + alpha);
      lp.set_b(current_row + 1, -bi_dot_p + alpha);

      current_row += 2;
    }

    //=========== Other LP parameters ===========
    lp.set_c(0, 1); // Minimize x[0]

    //=========== Solve =========================
    LP_solution solution = CGAL::solve_linear_program(lp, ET());
    bool ret = (solution.status() == CGAL::QP_OPTIMAL);

    return ret;
  }

  std::ostream &export_simplices_to_off(
    std::ostream & os, std::size_t &num_OFF_simplices,
    bool color_inconsistencies = false,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_red = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_green = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_blue = NULL)
    const
  {
    // If m_intrinsic_dim = 1, each point is output two times
    // (see export_vertices_to_off)
    num_OFF_simplices = 0;
    std::size_t num_maximal_simplices = 0;
    std::size_t num_inconsistent_maximal_simplices = 0;
    std::size_t num_inconsistent_stars = 0;
    typename Tr_container::const_iterator it_tr = m_triangulations.begin();
    typename Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for (std::size_t idx = 0 ; it_tr != it_tr_end ; ++it_tr, ++idx)
    {
      bool is_star_inconsistent = false;

      Triangulation const& tr    = it_tr->tr();
      Tr_vertex_handle center_vh = it_tr->center_vertex();

      if (&tr == NULL || tr.current_dimension() < m_intrinsic_dim)
        continue;

      // Color for this star
      std::stringstream color;
      //color << rand()%256 << " " << 100+rand()%156 << " " << 100+rand()%156;
      color << 128 << " " << 128 << " " << 128;

      // Gather the triangles here, with an int telling its color
      typedef std::vector<std::pair<Indexed_simplex, int> >
                                                          Star_using_triangles;
      Star_using_triangles star_using_triangles;

      // For each cell of the star
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        Indexed_simplex c = *it_inc_simplex;
        c.insert(idx);
        std::size_t num_vertices = c.size();
        ++num_maximal_simplices;

        int color_simplex = -1;// -1=no color, 0=yellow, 1=red, 2=green, 3=blue
        if (color_inconsistencies && !is_simplex_consistent(c))
        {
          ++num_inconsistent_maximal_simplices;
          color_simplex = 0;
          is_star_inconsistent = true;
        }
        else
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

        // If m_intrinsic_dim = 1, each point is output two times,
        // so we need to multiply each index by 2
        // And if only 2 vertices, add a third one (each vertex is duplicated in
        // the file when m_intrinsic dim = 2)
        if (m_intrinsic_dim == 1)
        {
          Indexed_simplex tmp_c;
          Indexed_simplex::iterator it = c.begin();
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
            Indexed_simplex triangle;
            Indexed_simplex::iterator it = c.begin();
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
        const Indexed_simplex &c = it_simplex->first;

        // Don't export infinite cells
        if (is_infinite(c))
          continue;
        
#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
        if (is_one_of_the_coord_far_from_origin(c, CGAL_TC_ALVAREZ_SURFACE_WINDOW, 2))
          continue;
#endif

        int color_simplex = it_simplex->second;

        std::stringstream sstr_c;

        Indexed_simplex::const_iterator it_point_idx = c.begin();
        for ( ; it_point_idx != c.end() ; ++it_point_idx)
        {
          sstr_c << *it_point_idx << " ";
        }

        // In order to have only one time each simplex, we only keep it
        // if the lowest index is the index of the center vertex
        // CJTODO: uncomment? but it only works if there's no inconsistencies
        /*if (*c.begin() != (m_intrinsic_dim == 1 ? 2*idx : idx)
            && color_simplex == -1)
          continue;*/

        os << 3 << " " << sstr_c.str();
        if (color_inconsistencies || p_simpl_to_color_in_red
            || p_simpl_to_color_in_green || p_simpl_to_color_in_blue)
        {
          switch (color_simplex)
          {
            case 0: os << " 255 255 0"; break;
            case 1: os << " 255 0 0"; break;
            case 2: os << " 0 255 0"; break;
            case 3: os << " 0 0 255"; break;
            default: os << " " << color.str(); break;
          }            
        }
        ++num_OFF_simplices;
        os << "\n";
      }
      if (is_star_inconsistent)
        ++num_inconsistent_stars;
    }

#ifdef CGAL_TC_VERBOSE
    std::cerr
      << "\n==========================================================\n"
      << "Export from list of stars to OFF:\n"
      << "  * Number of vertices: " << m_points.size() << "\n"
      << "  * Total number of maximal simplices: " << num_maximal_simplices 
      << "\n";
    if (color_inconsistencies)
    {
      std::cerr
        << "  * Number of inconsistent stars: "
        << num_inconsistent_stars << " ("
        << (m_points.size() > 0 ?
            100. * num_inconsistent_stars / m_points.size() : 0.) << "%)\n"
        << "  * Number of inconsistent maximal simplices: "
        << num_inconsistent_maximal_simplices << " ("
        << (num_maximal_simplices > 0 ?
            100. * num_inconsistent_maximal_simplices / num_maximal_simplices 
            : 0.) << "%)\n";
    }
    std::cerr << "==========================================================\n";
#endif

    return os;
  }

public:
  std::ostream &export_simplices_to_off(
    const Simplicial_complex &complex,
    std::ostream & os, std::size_t &num_OFF_simplices,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_red = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_green = NULL,
    std::set<Indexed_simplex > const *p_simpl_to_color_in_blue = NULL)
    const
  {
    typedef Simplicial_complex::Simplex                     Simplex;
    typedef Simplicial_complex::Simplex_set                 Simplex_set;

    // If m_intrinsic_dim = 1, each point is output two times
    // (see export_vertices_to_off)
    num_OFF_simplices = 0;
    std::size_t num_maximal_simplices = 0;

    typename Simplex_set::const_iterator it_s =
      complex.simplex_range().begin();
    typename Simplex_set::const_iterator it_s_end =
      complex.simplex_range().end();
    // For each simplex
    for ( ; it_s != it_s_end ; ++it_s)
    {
      Simplex c = *it_s;
      ++num_maximal_simplices;
      
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
      if (num_vertices < m_intrinsic_dim + 1)
        continue;

      // If m_intrinsic_dim = 1, each point is output two times,
      // so we need to multiply each index by 2
      // And if only 2 vertices, add a third one (each vertex is duplicated in
      // the file when m_intrinsic dim = 2)
      if (m_intrinsic_dim == 1)
      {
        Indexed_simplex tmp_c;
        Indexed_simplex::iterator it = c.begin();
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
          Indexed_simplex triangle;
          Indexed_simplex::iterator it = c.begin();
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
        if (is_infinite(*it_tri))
          continue;
        
#ifdef CGAL_TC_ALVAREZ_SURFACE_WINDOW
        if (is_one_of_the_coord_far_from_origin(*it_tri, CGAL_TC_ALVAREZ_SURFACE_WINDOW, 2))
          continue;
#endif

        os << 3 << " ";
        Indexed_simplex::const_iterator it_point_idx = it_tri->begin();
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

        ++num_OFF_simplices;
        os << "\n";
      }
    }

#ifdef CGAL_TC_VERBOSE
    std::cerr
      << "\n==========================================================\n"
      << "Export from complex to OFF:\n"
      << "  * Number of vertices: " << m_points.size() << "\n"
      << "  * Total number of maximal simplices: " << num_maximal_simplices 
      << "\n"
      << "==========================================================\n";
#endif

    return os;
  }

  // Return a pair<num_simplices, num_inconsistent_simplices>
  void check_correlation_between_inconsistencies_and_fatness() const
  {
    std::ofstream csv_consistent("output/correlation_consistent.csv"); // CJTODO DEBUG
    std::ofstream csv_inconsistent("output/correlation_inconsistent.csv"); // CJTODO DEBUG
    if (m_intrinsic_dim < 3)
    {
      std::cerr
        << "\n==========================================================\n"
        << "check_correlation_between_inconsistencies_and_fatness():\n"
        << "Intrinsic dimension should be >= 3.\n"
        << "==========================================================\n\n";
    }

    std::size_t num_consistent_simplices = 0;
    double sum_vol_edge_ratio_consistent = 0.;
    std::size_t num_inconsistent_simplices = 0;
    double sum_vol_edge_ratio_inconsistent = 0.;
    // For each triangulation
    for (std::size_t idx = 0 ; idx < m_points.size() ; ++idx)
    {
      // For each cell
      Star::const_iterator it_inc_simplex = m_stars[idx].begin();
      Star::const_iterator it_inc_simplex_end = m_stars[idx].end();
      for ( ; it_inc_simplex != it_inc_simplex_end ; ++it_inc_simplex)
      {
        // Don't check infinite cells
        if (is_infinite(*it_inc_simplex))
          continue;

        Indexed_simplex c = *it_inc_simplex;
        c.insert(idx); // Add the missing index

        double fatness = compute_simplex_fatness(c);
        
        if (!is_simplex_consistent(c))
        {
          ++num_inconsistent_simplices;
          sum_vol_edge_ratio_inconsistent += fatness;
          csv_inconsistent << fatness << "\n";
        }
        else
        {
          ++num_consistent_simplices;
          sum_vol_edge_ratio_consistent += fatness;
          csv_consistent << fatness << "\n";
        }
      }
    }

    double avg_vol_edge_ratio_inconsistent = 
      sum_vol_edge_ratio_inconsistent / num_inconsistent_simplices;
    double avg_vol_edge_ratio_consistent = 
      sum_vol_edge_ratio_consistent / num_consistent_simplices;
    
    std::cerr
      << "\n==========================================================\n"
      << "check_correlation_between_inconsistencies_and_fatness()\n"
      << "  * Avg. volume/longest_edge^d ratio of consistent simplices: " 
      << avg_vol_edge_ratio_consistent 
      << " (" << num_consistent_simplices << " simplices)\n"
      << "  * Avg. volume/longest_edge^d ratio of inconsistent simplices: " 
      << avg_vol_edge_ratio_inconsistent
      << " (" << num_inconsistent_simplices << " simplices)\n"
      << "==========================================================\n";
  }

private:
  const K              m_k;
  const int                 m_intrinsic_dim;
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
#ifdef CGAL_TC_EXPORT_NORMALS
  OS_container              m_orth_spaces;
#endif
  Tr_container              m_triangulations; // Contains the triangulations
                                              // and their center vertex
  Stars_container           m_stars;
  std::vector<FT>           m_squared_star_spheres_radii_incl_margin;
#ifdef CGAL_LINKED_WITH_TBB
  //std::vector<Tr_mutex>     m_tr_mutexes;
#endif

#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  Points                    m_points_for_tse;
  Points_ds                 m_points_ds_for_tse;
#endif

  mutable CGAL::Random      m_random_generator;

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
