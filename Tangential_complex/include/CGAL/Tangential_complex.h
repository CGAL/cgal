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

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Dimension.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Tangential_complex/utilities.h>
#include <CGAL/Tangential_complex/Point_cloud.h>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/point_generators_d.h>

#ifdef CGAL_TC_PROFILING
# include <CGAL/Mesh_3/Profiling_tools.h>
#endif

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

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/parallel_for.h>
//# include <tbb/queuing_mutex.h>
#endif

//#define CGAL_TC_EXPORT_NORMALS // Only for 3D surfaces (k=2, d=3)

namespace CGAL {

using namespace Tangential_complex_;

enum Fix_inconsistencies_status { TC_FIXED = 0, TIME_LIMIT_REACHED };

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
  typename Kernel,
  typename DimensionTag,
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
  typedef typename Kernel::Vector_d                   Vector;

  typedef Tr                                          Triangulation;
  typedef typename Triangulation::Geom_traits         Tr_traits;
  typedef typename Triangulation::Weighted_point      Tr_point;
  typedef typename Triangulation::Bare_point          Tr_bare_point;
  typedef typename Triangulation::Vertex_handle       Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle    Tr_full_cell_handle;
  typedef typename Tr_traits::Vector_d                Tr_vector;

  typedef std::vector<Vector>                         Tangent_space_basis;

  typedef std::vector<Point>                          Points;
  typedef std::vector<FT>                             Weights;
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

  typedef typename std::vector<Tangent_space_basis>   TS_container;
  typedef typename std::vector<Tr_and_VH>             Tr_container;
  typedef typename std::vector<Vector>                Vectors;
#ifdef CGAL_LINKED_WITH_TBB
  // CJTODO: test other mutexes
  // http://www.threadingbuildingblocks.org/docs/help/reference/synchronization/mutexes/mutex_concept.htm
  //typedef tbb::queuing_mutex                          Tr_mutex;
#endif

  // For transform_iterator
  static const Tr_point &vertex_handle_to_point(Tr_vertex_handle vh)
  {
    return vh->point();
  }

public:

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
    m_ambiant_dim(k.point_dimension_d_object()(*first)),
    m_points(first, last),
    m_points_ds(m_points)
#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
    , m_points_for_tse(first_for_tse, last_for_tse)
    , m_points_ds_for_tse(m_points_for_tse)
#endif
  {}

  /// Destructor
  ~Tangential_complex() {}

  std::size_t number_of_vertices()
  {
    return m_points.size();
  }

  void compute_tangential_complex()
  {
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif

    // We need to do that because we don't want the container to copy the
    // already-computed triangulations (while resizing) since it would
    // invalidate the vertex handles stored beside the triangulations
    m_triangulations.resize(m_points.size());
#ifdef CGAL_LINKED_WITH_TBB
    //m_tr_mutexes.resize(m_points.size());
#endif
    m_tangent_spaces.resize(m_points.size());
#ifdef CGAL_TC_PERTURB_WEIGHT
    m_weights.resize(m_points.size(), FT(0));
#endif
#ifdef CGAL_TC_PERTURB_POSITION
    m_translations.resize(m_points.size(),
                          m_k.construct_vector_d_object()(m_ambiant_dim));
#endif
#ifdef CGAL_TC_EXPORT_NORMALS
    m_normals.resize(m_points.size(),
                     m_k.construct_vector_d_object()(m_ambiant_dim));
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
    std::cerr << "Tangential complex computed in " << t.elapsed()
              << " seconds." << std::endl;
#endif
  }
  
  void estimate_intrinsic_dimension()
  {
    // Kernel functors
    typename Kernel::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();

    std::vector<FT> sum_eigen_values(m_ambiant_dim, FT(0));

    Points::const_iterator it_p = m_points.begin();
    Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for ( ; it_p != it_p_end ; ++it_p)
    {
      const Point &p = *it_p;
  
      KNS_range kns_range = m_points_ds.query_ANN(
        p, NUM_POINTS_FOR_PCA, false);

      //******************************* PCA *************************************

      // One row = one point
      Eigen::MatrixXd mat_points(NUM_POINTS_FOR_PCA, m_ambiant_dim);
      KNS_iterator nn_it = kns_range.begin();
      for (int j = 0 ;
            j < NUM_POINTS_FOR_PCA && nn_it != kns_range.end() ;
            ++j, ++nn_it)
      {
        for (int i = 0 ; i < m_ambiant_dim ; ++i)
          mat_points(j, i) = CGAL::to_double(coord(m_points[nn_it->first], i));
      }
      Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
      Eigen::MatrixXd cov = centered.adjoint() * centered;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

      // The eigenvectors are sorted in increasing order of their corresponding
      // eigenvalues
      Tangent_space_basis ts;
      for (int i = 0 ; i < m_ambiant_dim ; ++i)
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
        Compute_tangent_triangulation(*this, 
          true) //tangent_spaces_are_already_computed 
      );
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      for (std::size_t i = 0 ; i < m_points.size() ; ++i)
      {
        compute_tangent_triangulation(i, 
          true); // tangent_spaces_are_already_computed
      }
    }
    

#ifdef CGAL_TC_PROFILING
    std::cerr << "Tangential complex refreshed in " << t.elapsed()
              << " seconds." << std::endl;
#endif
  }

  // time_limit in seconds
  Fix_inconsistencies_status fix_inconsistencies(
    unsigned int &num_steps, 
    std::size_t &initial_num_inconsistent_local_tr,
    std::size_t &best_num_inconsistent_local_tr, 
    std::size_t &final_num_inconsistent_local_tr, 
    double time_limit = 0.)
  {
    Wall_clock_timer t;

    typename Kernel::Point_drop_weight_d drop_w = 
      m_k.point_drop_weight_d_object();

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
// CJTODO: the parallel version is not working for now
/*#ifdef CGAL_LINKED_WITH_TBB
      // Parallel
      if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
      {
        tbb::parallel_for(
          tbb::blocked_range<size_t>(0, m_triangulations.size()),
          Try_to_solve_inconsistencies_in_a_local_triangulation(*this)
        );
      }
      // Sequential
      else
#endif // CGAL_LINKED_WITH_TBB*/
      {
#ifdef CGAL_TC_PROFILING
        Wall_clock_timer t;
#endif
        for (std::size_t i = 0 ; i < m_triangulations.size() ; ++i)
        {
          num_inconsistent_local_tr += 
            (try_to_solve_inconsistencies_in_a_local_triangulation(i) ? 1 : 0);
        }
#ifdef CGAL_TC_PROFILING
    std::cerr << "Attempt to fix inconsistencies: " << t.elapsed()
              << " seconds." << std::endl;
#endif
      }

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
        << "  * BEFORE fix_inconsistencies:" << std::endl
        << "    - Total number of simplices in stars (incl. duplicates): "
        << stats_before.first << std::endl
        << "    - Num inconsistent simplices in stars (incl. duplicates): "
        << stats_before.second 
        << " (" << 100. * stats_before.second / stats_before.first << "%)"
        << std::endl
        << "  * Num inconsistent local triangulations: "
        << num_inconsistent_local_tr 
        << " (" << 100. * num_inconsistent_local_tr / m_points.size() << "%)" 
        << std::endl
        << std::endl
        << "  * AFTER fix_inconsistencies:" << std::endl
        << "    - Total number of simplices in stars (incl. duplicates): "
        << stats_after.first << std::endl
        << "    - Num inconsistent simplices in stars (incl. duplicates): "
        << stats_after.second << std::endl
        << "    - Percentage of inconsistencies: "
        << 100. * stats_after.second / stats_before.first << "%"
        << std::endl
        << "=========================================================="
        << std::endl;

      stats_before = stats_after;

#else // CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES
# ifdef CGAL_TC_VERBOSE
      std::cerr << std::endl
        << "=========================================================="
        << std::endl
        << "fix_inconsistencies():\n"
        << "  * " << m_points.size() << " vertices" << std::endl
        << "  * " << num_inconsistent_local_tr 
        << " (" << 100. * num_inconsistent_local_tr / m_points.size() << "%)"
        << " inconsistent triangulations encountered" << std::endl
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
      if (time_limit != 0 && t.elapsed() > time_limit)
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
    )
  {
    std::size_t num_simplices = 0;
    std::size_t num_inconsistent_simplices = 0;
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
        if (tr.is_infinite(*it_c)) // Don't check infinite cells
          continue;

        if (!is_simplex_consistent(*it_c, tr.current_dimension()))
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

  std::ostream &export_to_off(
    std::ostream & os,
    bool color_inconsistencies = false,
    std::set<std::set<std::size_t> > const* excluded_simplices = NULL,
    bool show_excluded_vertices_in_color = false)
  {
    if (m_points.empty())
      return os;

    const int ambient_dim = m_k.point_dimension_d_object()(*m_points.begin());
    if (ambient_dim < 2)
    {
      std::cerr << "Error: export_to_off => ambient dimension should be >= 2."
                << std::endl;
      os << "Error: export_to_off => ambient dimension should be >= 2."
         << std::endl;
      return os;
    }
    if (ambient_dim > 3)
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
    export_simplices_to_off(
      output, num_simplices, color_inconsistencies,
      excluded_simplices, show_excluded_vertices_in_color);

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

  bool check_if_all_simplices_are_in_the_ambient_delaunay(
    std::set<std::set<std::size_t> > * incorrect_simplices = NULL)
  {
    if (m_points.empty())
      return true;

    const int ambient_dim = m_k.point_dimension_d_object()(*m_points.begin());
    typedef Delaunay_triangulation<Kernel,
                                   Triangulation_data_structure<
                                     typename Kernel::Dimension,
                                     Triangulation_vertex<Kernel, Vertex_data>
                                   > >                        DT;
    typedef typename DT::Vertex_handle                        DT_VH;
    typedef typename DT::Finite_full_cell_const_iterator      FFC_it;
    typedef std::set<std::size_t>                             Indexed_simplex;

    //-------------------------------------------------------------------------
    // Build the ambient Delaunay triangulation
    // Then save its simplices into "amb_dt_simplices"
    //-------------------------------------------------------------------------

    DT ambient_dt(ambient_dim);
    for (std::size_t i=0; i<m_points.size(); ++i)
    {
      const Point& p = m_points[i];
      DT_VH vh = ambient_dt.insert(p);
      vh->data() = i;
    }

    std::set<Indexed_simplex> amb_dt_simplices;

    for (FFC_it cit = ambient_dt.finite_full_cells_begin() ;
         cit != ambient_dt.finite_full_cells_end() ; ++cit )
    {
      CGAL::Combination_enumerator<int> combi(
        m_intrinsic_dimension + 1, 0, ambient_dim + 1);

      for ( ; !combi.finished() ; ++combi)
      {
        Indexed_simplex simplex;
        for (int i = 0 ; i < m_intrinsic_dimension + 1 ; ++i)
          simplex.insert(cit.base()->vertex(combi[i])->data());

        amb_dt_simplices.insert(simplex);
      }
    }

    //-------------------------------------------------------------------------
    // Parse the TC and save its simplices into "stars_simplices"
    //-------------------------------------------------------------------------

    std::set<Indexed_simplex> stars_simplices;

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
        Indexed_simplex simplex;
        for (int i = 0 ; i < tr.current_dimension() + 1 ; ++i)
          simplex.insert((*it_c)->vertex(i)->data());

        stars_simplices.insert(simplex);
      }
    }

    //-------------------------------------------------------------------------
    // Check if simplices of "stars_simplices" are all in "amb_dt_simplices"
    //-------------------------------------------------------------------------

    std::set<Indexed_simplex> diff;
    if (!incorrect_simplices)
      incorrect_simplices = &diff;
    set_difference(stars_simplices.begin(),  stars_simplices.end(),
                   amb_dt_simplices.begin(), amb_dt_simplices.end(),
                   std::inserter(*incorrect_simplices,
                                 incorrect_simplices->begin()) );

    if (!incorrect_simplices->empty())
    {
      std::cerr
        << "ERROR check_if_all_simplices_are_in_the_ambient_delaunay:"
        << std::endl
        << "  Number of simplices in ambient DT: " << amb_dt_simplices.size()
        << std::endl
        << "  Number of unique simplices in TC stars: " << stars_simplices.size()
        << std::endl
        << "  Number of wrong simplices: " << incorrect_simplices->size()
        << std::endl;
      return false;
    }
    else
      return true;
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

  struct Tr_vertex_to_global_point
  {
    typedef Tr_vertex_handle argument_type;
    typedef Point result_type;

    Tr_vertex_to_global_point(Points const& points)
      : m_points(points) {}

    result_type operator()(argument_type const& vh) const
    {
      return m_points[vh->data()];
    }

  private:
    Points const& m_points;
  };

  struct Tr_vertex_to_bare_point
  {
    typedef Tr_vertex_handle argument_type;
    typedef Tr_bare_point result_type;

    Tr_vertex_to_bare_point(Tr_traits const& traits)
      : m_traits(traits) {}

    result_type operator()(argument_type const& vh) const
    {
      typename Tr_traits::Point_drop_weight_d pdw =
        m_traits.point_drop_weight_d_object();
      return pdw(vh->point());
    }

  private:
    Tr_traits const& m_traits;
  };

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for compute_tangential_complex function
  class Compute_tangent_triangulation
  {
    Tangential_complex & m_tc;
    bool m_tangent_spaces_are_already_computed;

  public:
    // Constructor
    Compute_tangent_triangulation(
      Tangential_complex &tc, bool tangent_spaces_are_already_computed = false)
    : m_tc(tc),
      m_tangent_spaces_are_already_computed(tangent_spaces_are_already_computed)
    {}

    // Constructor
    Compute_tangent_triangulation(const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
      {
        m_tc.compute_tangent_triangulation(
          i, m_tangent_spaces_are_already_computed);
      }
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  void compute_tangent_triangulation(
    std::size_t i, bool tangent_spaces_are_already_computed = false)
  {
    //std::cerr << "***********************************************" << std::endl;
    Triangulation &local_tr =
      m_triangulations[i].construct_triangulation(m_intrinsic_dimension);
    const Tr_traits &local_tr_traits = local_tr.geom_traits();
    Tr_vertex_handle &center_vertex = m_triangulations[i].center_vertex();

    // Kernel functor & objects
    typename Kernel::Difference_of_points_d k_diff_pts =
      m_k.difference_of_points_d_object();
    typename Kernel::Squared_distance_d k_sqdist =
      m_k.squared_distance_d_object();
    typename Kernel::Translated_point_d k_transl =
      m_k.translated_point_d_object();

    // Triangulation's traits functor & objects
    typename Tr_traits::Squared_distance_d sqdist =
      local_tr_traits.squared_distance_d_object();
    typename Tr_traits::Point_weight_d point_weight =
      local_tr_traits.point_weight_d_object();
    typename Tr_traits::Point_drop_weight_d drop_w =
      local_tr_traits.point_drop_weight_d_object();
    /*typename Tr_traits::Power_center_d power_center =
      local_tr_traits.power_center_d_object();*/ // CJTODO
    typename Get_functor<Tr_traits, Power_center_tag>::type power_center(local_tr_traits);

#ifdef CGAL_TC_PERTURB_POSITION
    const Point center_pt = k_transl(m_points[i], m_translations[i]);
#else
    const Point &center_pt = m_points[i];
#endif
    
    // Estimate the tangent space
    if (!tangent_spaces_are_already_computed)
    {
#ifdef CGAL_TC_EXPORT_NORMALS
      m_tangent_spaces[i] = compute_tangent_space(center_pt, &m_normals[i]);
#else
      m_tangent_spaces[i] = compute_tangent_space(center_pt);
#endif
    }

    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p)
    //***************************************************

    // Insert p
    Tr_point wp = local_tr_traits.construct_weighted_point_d_object()(
      local_tr_traits.construct_point_d_object()(m_intrinsic_dimension, ORIGIN),
#ifdef CGAL_TC_PERTURB_WEIGHT
      m_weights[i]
#else
      0
#endif
    );
    center_vertex = local_tr.insert(wp);
    center_vertex->data() = i;

    //const int NUM_NEIGHBORS = 150;
    //KNS_range ins_range = m_points_ds.query_ANN(center_pt, NUM_NEIGHBORS);
    INS_range ins_range = m_points_ds.query_incremental_ANN(center_pt);

    // While building the local triangulation, we keep the radius
    // of the sphere "star sphere" centered at "center_vertex"
    // and which contains all the
    // circumspheres of the star of "center_vertex"
    boost::optional<FT> star_sphere_squared_radius;

    // Insert points until we find a point which is outside "star shere"
    for (INS_iterator nn_it = ins_range.begin() ;
         nn_it != ins_range.end() ;
         ++nn_it)
    {
      std::size_t neighbor_point_idx = nn_it->first;

      // ith point = p, which is already inserted
      if (neighbor_point_idx != i)
      {
#ifdef CGAL_TC_PERTURB_POSITION
        const Point neighbor_pt = k_transl(
          m_points[neighbor_point_idx], m_translations[neighbor_point_idx]);
#else
        const Point &neighbor_pt = m_points[neighbor_point_idx];
#endif
#ifdef CGAL_TC_PERTURB_WEIGHT
        FT neighbor_weight = m_weights[neighbor_point_idx];
#else
        FT neighbor_weight = 0;
#endif

        // "4*m_sq_half_sparsity" because both points can be perturbed
        if (star_sphere_squared_radius
          && k_sqdist(center_pt, neighbor_pt)
             > *star_sphere_squared_radius + 4*m_sq_half_sparsity)
          break;

        Tr_point proj_pt = project_point_and_compute_weight(
          neighbor_pt, neighbor_weight, center_pt, m_tangent_spaces[i], 
          local_tr_traits);

        Tr_vertex_handle vh = local_tr.insert_if_in_star(proj_pt, center_vertex);
        //Tr_vertex_handle vh = local_tr.insert(proj_pt);
        if (vh != Tr_vertex_handle())
        {
          // CJTODO TEMP TEST
          /*if (star_sphere_squared_radius
            && k_sqdist(center_pt, neighbor_pt)
               > *star_sphere_squared_radius + 4*m_sq_half_sparsity)
            std::cout << "ARGGGGGGGH" << std::endl;*/

          vh->data() = neighbor_point_idx;

          // Let's recompute star_sphere_squared_radius
          if (local_tr.current_dimension() >= m_intrinsic_dimension)
          {
            star_sphere_squared_radius = 0;
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
                star_sphere_squared_radius = boost::none;
                break;
              }
              else
              {
                Tr_point c = power_center(
                  boost::make_transform_iterator(
                    cell->vertices_begin(), vertex_handle_to_point),
                  boost::make_transform_iterator(
                    cell->vertices_end(), vertex_handle_to_point));

                FT sq_power_sphere_diam = 4*point_weight(c);

                if (!star_sphere_squared_radius
                  || sq_power_sphere_diam > *star_sphere_squared_radius)
                {
                  star_sphere_squared_radius = sq_power_sphere_diam;
                }
              }
            }
          }
        }
        //std::cerr << star_sphere_squared_radius << std::endl;
      }
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

  Tangent_space_basis compute_tangent_space(const Point &p
#ifdef CGAL_TC_EXPORT_NORMALS
                                            , Vector *p_normal
#endif
                                            ) const
  {
    //******************************* PCA *************************************

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
    
#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
    KNS_range kns_range = m_points_ds_for_tse.query_ANN(
      p, NUM_POINTS_FOR_PCA, false);
    const Points &points_for_pca = m_points_for_tse;
#else
    KNS_range kns_range = m_points_ds.query_ANN(p, NUM_POINTS_FOR_PCA, false);
    const Points &points_for_pca = m_points;
#endif

    // One row = one point
    Eigen::MatrixXd mat_points(NUM_POINTS_FOR_PCA, m_ambiant_dim);
    KNS_iterator nn_it = kns_range.begin();
    for (int j = 0 ;
         j < NUM_POINTS_FOR_PCA && nn_it != kns_range.end() ;
         ++j, ++nn_it)
    {
      for (int i = 0 ; i < m_ambiant_dim ; ++i)
        mat_points(j, i) = CGAL::to_double(coord(points_for_pca[nn_it->first], i));
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    Tangent_space_basis ts;
    for (int i = m_ambiant_dim - 1 ; 
         i >= m_ambiant_dim - m_intrinsic_dimension ; 
         --i)
    {
      ts.push_back(constr_vec(
        m_ambiant_dim,
        eig.eigenvectors().col(i).data(),
        eig.eigenvectors().col(i).data() + m_ambiant_dim));
    }
#ifdef CGAL_TC_EXPORT_NORMALS
    *p_normal = constr_vec(
        m_ambiant_dim,
        eig.eigenvectors().col(m_ambiant_dim - m_intrinsic_dimension - 1).data(),
        eig.eigenvectors().col(m_ambiant_dim - m_intrinsic_dimension - 1).data() + m_ambiant_dim);
#endif

    //*************************************************************************

    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, FT(1)/sqrt(sqlen(n)));
    //std::cerr << "IP = " << inner_pdct(n, ts[0]) << " & " << inner_pdct(n, ts[1]) << std::endl;

    return compute_gram_schmidt_basis(ts, m_k);

    
    // CJTODO: this is only for a sphere in R^3
    /*double tt1[3] = {-p[1] - p[2], p[0], p[0]};
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

    return ts;*/

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
    //return compute_gram_schmidt_basis(ts, m_k);
    */
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_bare_point project_point(const Point &p, const Point &origin,
                         const Tangent_space_basis &ts) const
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
      // Compute the inner product p * ts[i]
      Vector v = diff_points(p, origin);
      FT coord = inner_pdct(v, ts[i]);
      coords.push_back(coord);
    }

    return Tr_bare_point(m_intrinsic_dimension, coords.begin(), coords.end());
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_point project_point_and_compute_weight(
    const Point &p, FT w, const Point &origin, const Tangent_space_basis &ts,
    const Tr_traits &tr_traits) const
  {
    const int point_dim = m_k.point_dimension_d_object()(p);
    typename Kernel::Scalar_product_d inner_pdct =
      m_k.scalar_product_d_object();
    typename Kernel::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();
    typename Kernel::Construct_cartesian_const_iterator_d ccci =
      m_k.construct_cartesian_const_iterator_d_object();

    Vector v = diff_points(p, origin);

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(ccci(origin), ccci(origin, 0));
    coords.reserve(m_intrinsic_dimension);
    for (std::size_t i = 0 ; i < m_intrinsic_dimension ; ++i)
    {
      // Compute the inner product p * ts[i]
      FT coord = inner_pdct(v, ts[i]);
      coords.push_back(coord);

      // p_proj += coord * v;
      for (int j = 0 ; j < point_dim ; ++j)
        p_proj[j] += coord * ts[i][j];
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
  bool is_simplex_consistent(Tr_full_cell_handle fch, int cur_dim)
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
  bool is_simplex_consistent(std::set<std::size_t> const& simplex)
  {
    int cur_dim_plus_1 = static_cast<int>(simplex.size());

    // Check if the simplex is in the stars of all its vertices
    std::set<std::size_t>::const_iterator it_point_idx = simplex.begin();
    // For each point p of the simplex, we parse the incidents cells of p
    // and we check if "simplex" is among them
    for ( ; it_point_idx != simplex.end() ; ++it_point_idx)
    {
      std::size_t point_idx = *it_point_idx;
      if (point_idx == std::numeric_limits<std::size_t>::max())
        continue;
      Triangulation const& tr = m_triangulations[point_idx].tr();
      Tr_vertex_handle center_vh = m_triangulations[point_idx].center_vertex();

      std::vector<Tr_full_cell_handle> incident_cells;
      tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

      typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
      typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end= incident_cells.end();
      // For each cell
      bool found = false;
      for ( ; !found && it_c != it_c_end ; ++it_c)
      {
        std::set<std::size_t> cell;
        for (int i = 0 ; i < cur_dim_plus_1 ; ++i)
          cell.insert((*it_c)->vertex(i)->data());
        if (cell == simplex)
          found = true;
      }

      if (!found)
        return false;
    }

    return true;
  }

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for try_to_solve_inconsistencies_in_a_local_triangulation function
  class Try_to_solve_inconsistencies_in_a_local_triangulation
  {
    Tangential_complex & m_tc;

  public:
    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
      Tangential_complex &tc)
    : m_tc(tc)
    {}

    // Constructor
    Try_to_solve_inconsistencies_in_a_local_triangulation(
      const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
        m_tc.try_to_solve_inconsistencies_in_a_local_triangulation(i);
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  void perturb(std::size_t point_idx)
  {
    CGAL::Random rng;

    // Perturb the weight?
#ifdef CGAL_TC_PERTURB_WEIGHT
    m_weights[point_idx] = rng.get_double(0., m_sq_half_sparsity);
#endif

    // Perturb the position?
#ifdef CGAL_TC_PERTURB_POSITION
# ifdef CGAL_TC_PERTURB_POSITION_GLOBAL
    typename Kernel::Point_to_vector_d k_pt_to_vec =
      m_k.point_to_vector_d_object();
    CGAL::Random_points_on_sphere_d<Point> 
      tr_point_on_sphere_generator(m_ambiant_dim, 1);

    m_translations[point_idx] = k_scaled_vec(k_pt_to_vec(
      *tr_point_on_sphere_generator++), m_half_sparsity);
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
      tr_point_on_sphere_generator(m_intrinsic_dimension, 1);

    Tr_point local_random_transl =
      local_tr_traits.construct_weighted_point_d_object()(
        *tr_point_on_sphere_generator++, 0);
    Vector &global_transl = m_translations[point_idx];
    global_transl = k_constr_vec(m_ambiant_dim);
    const Tangent_space_basis &tsb = m_tangent_spaces[point_idx];
    for (int i = 0 ; i < m_intrinsic_dimension ; ++i)
    {
      global_transl = k_transl(
        global_transl, 
        k_scaled_vec(tsb[i], m_half_sparsity*coord(local_random_transl, i))
      );
    }
# endif
#endif // CGAL_TC_PERTURB_POSITION
  }

  bool try_to_solve_inconsistencies_in_a_local_triangulation(
                                                          std::size_t tr_index)
  {
    bool is_inconsistent = false;

#ifdef CGAL_LINKED_WITH_TBB
    //Tr_mutex::scoped_lock lock(m_tr_mutexes[tr_index]);
#endif
    
    Triangulation const& tr    = m_triangulations[tr_index].tr();
    Tr_vertex_handle center_vh = m_triangulations[tr_index].center_vertex();
    const Tr_traits &local_tr_traits = tr.geom_traits();
    int cur_dim = tr.current_dimension();

    std::vector<Tr_full_cell_handle> incident_cells;
    tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

    typename std::vector<Tr_full_cell_handle>::const_iterator it_c =
                                                        incident_cells.begin();
    typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end=
                                                        incident_cells.end();
    // For each cell
    for ( ; it_c != it_c_end ; ++it_c)
    {
      if (tr.is_infinite(*it_c)) // Don't check infinite cells
        continue;

//*****************************************************************************
// STRATEGY 1: perturb all the points of the first inconsistent simplex
//*****************************************************************************
#ifdef CGAL_TC_PERTURB_THE_SIMPLEX_ONLY

      std::set<std::size_t> c;
      for (int i = 0 ; i < tr.current_dimension() + 1 ; ++i)
      {
        std::size_t data = (*it_c)->vertex(i)->data();
        c.insert(data);
      }

      // Inconsistent?
      if (!is_simplex_consistent(c))
      {
        is_inconsistent = true;

        for (std::set<std::size_t>::iterator it=c.begin(); it!=c.end(); ++it)
          perturb(*it);
        
# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time (incident_cells is not
        // valid anymore here)
        break;
      }

//*****************************************************************************
// STRATEGY 2: perturb the center point only
//*****************************************************************************
#elif defined(CGAL_TC_PERTURB_THE_CENTER_VERTEX_ONLY)
      if (!is_simplex_consistent(*it_c, cur_dim))
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

        // We will try the other cells next time (incident_cells is not
        // valid anymore here)
        break;
      }
      
//*****************************************************************************
// STRATEGY 3: perturb all the points of the 1-star
//*****************************************************************************
#elif defined(CGAL_TC_PERTURB_THE_1_STAR)

      // Inconsistent?
      if (!is_simplex_consistent(*it_c, cur_dim))
      {
        is_inconsistent = true;
        
        std::set<std::size_t> c;
      
        typename std::vector<Tr_full_cell_handle>::const_iterator it_c2 =
                                                        incident_cells.begin();
        // For each cell
        for ( ; it_c2 != it_c_end ; ++it_c2)
        {
          for (int i = 0 ; i < tr.current_dimension() + 1 ; ++i)
          {
            std::size_t data = (*it_c2)->vertex(i)->data();
            c.insert(data);
          }
        }
        
        for (std::set<std::size_t>::iterator it=c.begin(); it!=c.end(); ++it)
          perturb(*it);
        
# if !defined(CGAL_TC_GLOBAL_REFRESH)
        refresh_tangential_complex();
# endif

        // We will try the other cells next time (incident_cells is not
        // valid anymore here)
        break;
      }

//*****************************************************************************
// STRATEGY 4: perturb the k + 1 + CGAL_TC_NUMBER_OF_ADDITIONNAL_PERTURBED_POINTS
// closest points (to the power center of first the inconsistent cell)
//*****************************************************************************
#else
      // Inconsistent?
      if (!is_simplex_consistent(*it_c, cur_dim))
      {
        is_inconsistent = true;

        // Get the k + 1 + CGAL_TC_NUMBER_OF_ADDITIONNAL_PERTURBED_POINTS
        // closest points

        std::vector<Tr_point> simplex_pts;
        for (int i = 0 ; i < cur_dim + 1 ; ++i)
          simplex_pts.push_back((*it_c)->vertex(i)->point());

        //typename Tr_traits::Power_center_d power_center =
        //  local_tr_traits.power_center_d_object(); // CJTODO
        typename Get_functor<Tr_traits, Power_center_tag>::type power_center(local_tr_traits);
        typename Tr_traits::Compute_coordinate_d coord =
          local_tr_traits.compute_coordinate_d_object();
        
        typename Kernel::Translated_point_d k_transl =
          m_k.translated_point_d_object();
        typename Kernel::Scaled_vector_d k_scaled_vec =
          m_k.scaled_vector_d_object();


        Tr_point local_center = power_center(simplex_pts.begin(), simplex_pts.end());
        Point global_center = m_points[tr_index];
        const Tangent_space_basis &tsb = m_tangent_spaces[tr_index];
        for (int i = 0 ; i < m_intrinsic_dimension ; ++i)
        {
          global_center = k_transl(
            global_center, 
            k_scaled_vec(tsb[i], coord(local_center, i)));
        }
        
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

        // We will try the other cells next time (incident_cells is not
        // valid anymore here)
        break;
      }

#endif // CGAL_TC_PERTURB_THE_SIMPLEX_ONLY
    }

    return is_inconsistent;
  }

  std::ostream &export_vertices_to_off(
    std::ostream & os, std::size_t &num_vertices)
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
    const int ambient_dim = m_k.point_dimension_d_object()(*m_points.begin());

    // Kernel functors
    typename Kernel::Compute_coordinate_d coord =
      m_k.compute_coordinate_d_object();

    int num_coords = min(ambient_dim, 3);
#ifdef CGAL_TC_EXPORT_NORMALS
    Vectors::const_iterator it_n = m_normals.begin();
#endif
    typename Points::const_iterator it_p = m_points.begin();
    typename Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for ( ; it_p != it_p_end ; ++it_p)
    {
      for (int ii = 0 ; ii < N ; ++ii)
      {
        int i = 0;
        for ( ; i < num_coords ; ++i)
          os << CGAL::to_double(coord(*it_p, i)) << " ";
        if (i == 2)
          os << "0";

#ifdef CGAL_TC_EXPORT_NORMALS
        for (i = 0 ; i < num_coords ; ++i)
          os << " " << CGAL::to_double(coord(*it_n, i));
#endif
        os << std::endl;
      }
#ifdef CGAL_TC_EXPORT_NORMALS
      ++it_n;
#endif
    }

    num_vertices = N*m_points.size();
    return os;
  }

  std::ostream &export_simplices_to_off(
    std::ostream & os, std::size_t &num_simplices,
    bool color_inconsistencies = false,
    std::set<std::set<std::size_t> > const* excluded_simplices = NULL,
    bool show_excluded_vertices_in_color = false)
  {
    // If m_intrinsic_dimension = 1, each point is output two times
    // (see export_vertices_to_off)
    int factor = (m_intrinsic_dimension == 1 ? 2 : 1);
    int OFF_simplices_dim =
      (m_intrinsic_dimension == 1 ? 3 : m_intrinsic_dimension + 1);
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

      if (tr.current_dimension() < m_intrinsic_dimension)
        continue;

      // Color for this star
      std::stringstream color;
      //color << rand()%256 << " " << 100+rand()%156 << " " << 100+rand()%156;
      color << 128 << " " << 128 << " " << 128;

      std::vector<Tr_full_cell_handle> incident_cells;
      tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

      typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
      typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end= incident_cells.end();
      // For each cell
      for ( ; it_c != it_c_end ; ++it_c)
      {
        if (tr.is_infinite(*it_c)) // Don't export infinite cells
          continue;

        if (color_inconsistencies || excluded_simplices)
        {
          std::set<std::size_t> c;
          std::stringstream sstr_c;
          std::size_t data;
          for (int i = 0 ; i < m_intrinsic_dimension + 1 ; ++i)
          {
            data = (*it_c)->vertex(i)->data();
            sstr_c << data*factor << " ";
            c.insert(data);
          }

          bool is_simpl_consistent = is_simplex_consistent(c);
          if (!is_simpl_consistent) 
            is_star_inconsistent = true;

          // See export_vertices_to_off
          if (m_intrinsic_dimension == 1)
            sstr_c << (data*factor + 1) << " ";

          // In order to have only one time each simplex, we only keep it
          // if the lowest index is the index of the center vertex
          if (*c.begin() != idx && is_simpl_consistent)
            continue;

          bool excluded =
            (excluded_simplices
            && excluded_simplices->find(c) != excluded_simplices->end());

          if (!excluded)
          {
            os << OFF_simplices_dim << " " << sstr_c.str() << " ";
            if (color_inconsistencies && is_simpl_consistent)
              os << color.str();
            else
            {
              os << "255 0 0";
              ++num_inconsistent_simplices;
            }
            ++num_simplices;
          }
          else if (show_excluded_vertices_in_color)
          {
            os << OFF_simplices_dim << " "
               << sstr_c.str() << " "
               << "0 0 255";
            ++num_simplices;
          }
        }
        else
        {
          os << OFF_simplices_dim << " ";
          std::size_t min_index = std::numeric_limits<std::size_t>::max();
          std::size_t data;
          std::stringstream sstr_c;
          for (int i = 0 ; i < m_intrinsic_dimension + 1 ; ++i)
          {
            data = (*it_c)->vertex(i)->data();
            if (data < min_index)
              min_index = data;
            sstr_c << data*factor << " ";
          }
          
          // In order to have only one time each simplex, we only keep it
          // if the lowest index is the index of the center vertex
          if (min_index != idx)
            continue;

          // See export_vertices_to_off
          if (m_intrinsic_dimension == 1)
            sstr_c << (data*factor + 1) << " ";

          os << sstr_c.str();
          ++num_simplices;
        }

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

private:
  const Kernel              m_k;
  const int                 m_intrinsic_dimension;
  const double              m_half_sparsity;
  const double              m_sq_half_sparsity;
  const int                 m_ambiant_dim;
  
  Points                    m_points;
#ifdef CGAL_TC_PERTURB_WEIGHT
  Weights                   m_weights;
#endif
#ifdef CGAL_TC_PERTURB_POSITION
  Vectors                   m_translations;
#endif

  Points_ds                 m_points_ds;
  TS_container              m_tangent_spaces;
  Tr_container              m_triangulations; // Contains the triangulations
                                              // and their center vertex
#ifdef CGAL_LINKED_WITH_TBB
  //std::vector<Tr_mutex>     m_tr_mutexes;
#endif
#ifdef CGAL_TC_EXPORT_NORMALS
  Vectors                   m_normals;
#endif

#ifdef USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  Points                    m_points_for_tse;
  Points_ds                 m_points_ds_for_tse;
#endif

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
