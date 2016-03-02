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


#ifndef MESH_D_H
#define MESH_D_H

#include <CGAL/Tangential_complex/config.h>
#include <CGAL/Tangential_complex/Simplicial_complex.h>
#include <CGAL/Tangential_complex/utilities.h>
#include <CGAL/Tangential_complex/Point_cloud.h>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Dimension.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/point_generators_d.h>
# include <CGAL/Mesh_3/Profiling_tools.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <iterator>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/parallel_for.h>
#endif

// choose exact integral type for QP solver
// (Gmpzf is not thread-safe)
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
//#define CGAL_QP_NO_ASSERTIONS // CJTODO: NECESSARY? http://doc.cgal.org/latest/QP_solver/group__PkgQPSolverFunctions.html#ga1fefbd0436aca0e281f88e8e6cd8eb74

namespace CGAL {

using namespace Tangential_complex_;

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

/// The class Mesh_d represents a d-dimensional mesh
template <
  typename Kernel_, // ambiant kernel
  typename Local_kernel_, // local kernel (intrinsic dimension)
  typename Concurrency_tag = CGAL::Parallel_tag,
  typename Tr = Delaunay_triangulation
  <
    Kernel_,
    Triangulation_data_structure
    <
      typename Kernel_::Dimension,
      Triangulation_vertex<Kernel_, Vertex_data>,
      Triangulation_full_cell<Kernel_>
    >
  >
>
class Mesh_d
{
  typedef Kernel_                                     K;
  typedef typename K::FT                              FT;
  typedef typename K::Point_d                         Point;
  typedef typename K::Weighted_point_d                Weighted_point;
  typedef typename K::Vector_d                        Vector;

  typedef Local_kernel_                               LK;
  typedef typename LK::Point_d                        Local_point;
  typedef typename LK::Weighted_point_d               Local_weighted_point;

  typedef Tr                                          Triangulation;
  typedef typename Triangulation::Vertex_handle       Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle    Tr_full_cell_handle;
  typedef typename Triangulation::Finite_full_cell_const_iterator 
                                            Tr_finite_full_cell_const_iterator;

  typedef std::vector<Point>                          Points;

  typedef Point_cloud_data_structure<K, Points>       Points_ds;
  typedef typename Points_ds::KNS_range               KNS_range;
  typedef typename Points_ds::KNS_iterator            KNS_iterator;
  typedef typename Points_ds::INS_range               INS_range;
  typedef typename Points_ds::INS_iterator            INS_iterator;
  
  typedef std::set<Tr_vertex_handle>                  Vertex_set;
  typedef std::set<std::size_t>                       Indexed_simplex;

public:
  typedef Basis<K>                                    Tangent_space_basis;
  typedef Basis<K>                                    Orthogonal_space_basis;
  typedef std::vector<Tangent_space_basis>            TS_container;
  typedef std::vector<Orthogonal_space_basis>         OS_container;

private:

  // For transform_iterator
  static const Point &vertex_handle_to_point(Tr_vertex_handle vh)
  {
    return vh->point();
  }
  // For transform_iterator
  static std::size_t vertex_handle_to_index(Tr_vertex_handle vh)
  {
    return vh->data();
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
  Mesh_d(InputIterator first, InputIterator last,
         double sparsity, int intrinsic_dimension,
#ifdef CGAL_MESH_D_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
         InputIterator first_for_tse, InputIterator last_for_tse,
#endif
         const K &k = K(),
         const LK &lk = LK()
         )
  : m_k(k)
  , m_lk(lk)
  , m_intrinsic_dim(intrinsic_dimension)
  , m_half_sparsity(0.5*sparsity)
  , m_sq_half_sparsity(m_half_sparsity*m_half_sparsity)
  , m_ambient_dim(k.point_dimension_d_object()(*first))
  , m_points(first, last)
  , m_points_ds(m_points)
  , m_are_tangent_spaces_computed(m_points.size(), false)
  , m_tangent_spaces(m_points.size(), Tangent_space_basis())
  , m_orth_spaces(m_points.size(), Orthogonal_space_basis())
#ifdef CGAL_MESH_D_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  , m_points_for_tse(first_for_tse, last_for_tse)
  , m_points_ds_for_tse(m_points_for_tse)
#endif
  {
    if (sparsity <= 0.)
      std::cerr << "!Warning! Sparsity should be > 0\n";
  }

  /// Destructor
  ~Mesh_d()
  {
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

  Simplicial_complex const& complex()
  {
    return m_complex;
  }

  void set_tangent_planes(
    const TS_container& tangent_spaces, const OS_container& orthogonal_spaces)
  {
    CGAL_assertion(m_points.size() == tangent_spaces.size()
                   && m_points.size() == orthogonal_spaces.size());
    m_tangent_spaces = tangent_spaces;
    m_orth_spaces = orthogonal_spaces;
    for(std::size_t i=0; i<m_points.size(); ++i)
      m_are_tangent_spaces_computed[i] = true;
  }

  void compute_mesh(std::vector<Indexed_simplex> *p_uncertain_simplices = NULL)
  {
#if defined(CGAL_MESH_D_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    Wall_clock_timer t;
#endif

    typedef CGAL::Epick_d<Dynamic_dimension_tag>        Orth_K;
    typedef typename Orth_K::Point_d                    Orth_K_point;
    typedef typename Orth_K::Vector_d                   Orth_K_vector;

    int orth_dim = m_ambient_dim - m_intrinsic_dim;
    Orth_K orth_k(orth_dim);

    Get_functor<Orth_K, Sum_of_vectors_tag>::type orth_k_sum_vecs(orth_k);

    // Tangent and orthogonal spaces
    for (std::size_t i = 0; i < m_points.size(); ++i)
    {
      if (!m_are_tangent_spaces_computed[i])
      {
        m_tangent_spaces[i] = compute_tangent_space(
          m_points[i], i, true/*normalize*/, &m_orth_spaces[i]);
      }
    }
#if defined(CGAL_MESH_D_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    std::cerr << "Tangent & orthogonal spaces computed in " << t.elapsed() 
      << " seconds." << std::endl;
    t.reset();
#endif

    // Ambiant DT
    Tr tr(m_ambient_dim);
    for (std::size_t i = 0; i < m_points.size(); ++i)
    {
      Tr_vertex_handle vh = tr.insert(m_points[i]);
      vh->data() = i;
    }
#if defined(CGAL_MESH_D_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    std::cerr << "Ambient DT computed in " << t.elapsed() 
      << " seconds." << std::endl;
    t.reset();
#endif

    // Extract complex
    // CJTODO: avoid duplicates + parallelize
    /*m_complex.clear();
    for (Tr_finite_full_cell_const_iterator cit = tr.finite_full_cells_begin() ;
      cit != tr.finite_full_cells_end() ; ++cit)
    {
      // Enumerate k-dim simplices
      CGAL::Combination_enumerator<int> combi(
        m_intrinsic_dim + 1, 0, m_ambient_dim + 1);

      for (; !combi.finished() ; ++combi)
      {
        Indexed_simplex simplex;
        std::vector<Tr_vertex_handle> vertices;
        for (int i = 0 ; i < m_intrinsic_dim + 1 ; ++i)
        {
          vertices.push_back(cit->vertex(combi[i]));
          simplex.insert(cit->vertex(combi[i])->data());
        }

        Vertex_set Q_set =
          get_common_neighbor_vertices(tr, vertices);
        // Convert it to a vector, because otherwise, the result of
        // boost::adaptors::transform does not provide size()
        std::vector<Tr_vertex_handle> Q(Q_set.begin(), Q_set.end());
        bool intersect = does_voronoi_face_and_tangent_subspace_intersect(
          boost::adaptors::transform(vertices, vertex_handle_to_point),
          boost::adaptors::transform(Q, vertex_handle_to_point),
          m_orth_spaces[*simplex.begin()]); // CJTODO: remplacer simplex[0] par un truc plus intelligent
        if (intersect)
          m_complex.add_simplex(simplex, false);
      }
    }*/

    Get_functor<K, Sum_of_vectors_tag>::type sum_vecs(m_k);

    // Extract complex
    // CJTODO: parallelize
    m_complex.clear();
    typedef std::map<Vertex_set, Vertex_set> K_faces_and_neighbor_vertices;
    K_faces_and_neighbor_vertices k_faces_and_neighbor_vertices =
      get_k_faces_and_neighbor_vertices(tr);
    m_complex.clear();
    for (K_faces_and_neighbor_vertices::const_iterator k_face_and_nghb_it = 
      k_faces_and_neighbor_vertices.begin() ;
      k_face_and_nghb_it != k_faces_and_neighbor_vertices.end() ; 
      ++k_face_and_nghb_it)
    {
      Vertex_set const& k_face = k_face_and_nghb_it->first;
      Vertex_set const& neighbor_vertices = k_face_and_nghb_it->second;
      // Convert it to a vector, because otherwise, the result of
      // boost::adaptors::transform does not provide size()
      std::vector<Tr_vertex_handle> kf(k_face.begin(), k_face.end());
      std::vector<Tr_vertex_handle> nghb(
        neighbor_vertices.begin(), neighbor_vertices.end());

      bool keep_it = false;
      bool is_uncertain = false;
#ifdef CGAL_MESH_D_FILTER_BY_TESTING_ALL_VERTICES_TANGENT_PLANES
      int num_intersections = 0;
      for (auto vh : kf) // CJTODO C++11
      {
        bool intersect = does_voronoi_face_and_tangent_subspace_intersect(
          boost::adaptors::transform(kf, vertex_handle_to_point),
          boost::adaptors::transform(nghb, vertex_handle_to_point),
          m_orth_spaces[vh->data()]);
        if (intersect)
          ++num_intersections;
      }

      if (num_intersections >= 1)
      {
        keep_it = true;
        if (num_intersections < m_intrinsic_dim + 1)
          is_uncertain = true;
      }
#else
      // Compute the intersection with all tangent planes of the vertices
      // of the k-face then compute a weighted barycenter
      FT sum_weights = 0;
      Vector weighted_sum_of_inters =
        orth_k.construct_vector_d_object()(m_ambient_dim);
      for (auto vh : kf) // CJTODO C++11
      {
        Point intersection = 
          compute_aff_of_voronoi_face_and_tangent_subspace_intersection(
            boost::adaptors::transform(kf, vertex_handle_to_point),
#ifdef CGAL_MESH_D_USE_LINEAR_PROG_TO_COMPUTE_INTERSECTION
            m_orth_spaces[vh->data()]);
#else
            m_tangent_spaces[vh->data()]);
#endif

        FT weight = 
          FT(1) / m_k.squared_distance_d_object()(vh->point(), intersection);
        sum_weights += weight;
        weighted_sum_of_inters = sum_vecs(
          weighted_sum_of_inters,
          m_k.scaled_vector_d_object()(
            m_k.point_to_vector_d_object()(intersection), weight));
      }

      if (sum_weights > 0)
      {
        // Compute the weighted barycenter
        Point avg_inters = m_k.vector_to_point_d_object()(
          m_k.scaled_vector_d_object()(
            weighted_sum_of_inters, FT(1) / sum_weights));

        //****************************************************************
        // Translate the tangent plane to be closer to the surface
        //****************************************************************

        // Estimate the normal subspace at point "avg_inters"
        Tangent_space_basis tsb;
        Orthogonal_space_basis osb;
        tsb = compute_tangent_space(
          avg_inters, std::numeric_limits<std::size_t>::max(), 
          true/*normalize*/, &osb);

        unsigned int num_points_for_nghb_query = static_cast<unsigned int>(
          std::pow(BASE_VALUE_FOR_PCA, m_intrinsic_dim));
        KNS_range kns_range = m_points_ds.query_ANN(
          avg_inters, num_points_for_nghb_query, false);

        FT sum_weights_of_barycenter = 0;
        Orth_K_vector weighted_sum_of_neighbors =
          orth_k.construct_vector_d_object()(orth_dim);
        KNS_iterator nn_it = kns_range.begin();
        for (unsigned int j = 0 ;
          j < num_points_for_nghb_query && nn_it != kns_range.end() ;
          ++j, ++nn_it)
        {
          Point const& nghb = m_points[nn_it->first];

          Orth_K_point proj_nghb = project_point(nghb, osb, orth_k, &avg_inters);

          FT weight =
            FT(1) / m_k.squared_distance_d_object()(nghb, avg_inters);
          sum_weights_of_barycenter += weight;
          weighted_sum_of_neighbors = orth_k_sum_vecs(
            weighted_sum_of_neighbors,
            orth_k.scaled_vector_d_object()(
              orth_k.point_to_vector_d_object()(proj_nghb), weight));
        }

        //sum_weights_of_barycenter *= 8; // CJTODO TEMP

        // Compute the weighted barycenter
        Point translated_origin = unproject_point(
          orth_k.vector_to_point_d_object()(
            orth_k.scaled_vector_d_object()(
              weighted_sum_of_neighbors, FT(1) / sum_weights_of_barycenter)),
          osb,
          orth_k,
          &avg_inters);

        Point translated_inters =
          compute_aff_of_voronoi_face_and_tangent_subspace_intersection(
            boost::adaptors::transform(kf, vertex_handle_to_point),
            tsb, &translated_origin);

        // CJTODO TEMP: exact intersection for the sphere of radius 1
        //translated_inters = orth_k.vector_to_point_d_object()(
        //  normalize_vector(orth_k.point_to_vector_d_object()(avg_inters), m_k));

        //****************************************************************
        // Keep the simplex or not?
        //****************************************************************
        
        keep_it = true;
        // Check if the averaged intersection "i" is inside the Voronoi cell, i.e.
        // for each common neighbor qi, (i - p0)² <= (i - qi)²
        Point const& p0 = (*k_face.begin())->point();
        for (auto neighbor_vh : nghb) // CJTODO: C++11
        {
          // If a neighbor is closer than p0 to the average intersection, then
          // the average intersection is inside another Voronoi cell
          if (m_k.squared_distance_d_object()(p0, translated_inters) >
            m_k.squared_distance_d_object()(neighbor_vh->point(), translated_inters))
          {
            keep_it = false;
            break;
          }
        }
      }
#endif

      if (keep_it)
      {
        Indexed_simplex s(
          boost::make_transform_iterator(kf.begin(), vertex_handle_to_index),
          boost::make_transform_iterator(kf.end(), vertex_handle_to_index));
        m_complex.add_simplex(s, false);

        if (is_uncertain && p_uncertain_simplices)
          p_uncertain_simplices->push_back(s);
      }
    }

#if defined(CGAL_MESH_D_PROFILING) && defined(CGAL_LINKED_WITH_TBB)
    std::cerr << "Mesh extracted in computed in " << t.elapsed()
      << " seconds." << std::endl;
    t.reset();
#endif
  }

  void display_stats()
  {
    m_complex.display_stats();
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

  bool is_infinite(Indexed_simplex const& s) const
  {
    return *s.rbegin() == std::numeric_limits<std::size_t>::max();
  }

  // If i = std::numeric_limits<std::size_t>::max(), it is ignored
  Tangent_space_basis compute_tangent_space(
      const Point &p
    , const std::size_t i = std::numeric_limits<std::size_t>::max()
    , bool normalize_basis = true
    , Orthogonal_space_basis *p_orth_space_basis = NULL)
  {
#ifdef CGAL_MESH_D_COMPUTE_TANGENT_PLANES_FOR_SPHERE_2

    double tt[2] = {p[1], -p[0]};
    Vector t(2, &tt[0], &tt[2]);

    // Normalize t1 and t2
    typename K::Squared_length_d sqlen      = m_k.squared_length_d_object();
    typename K::Scaled_vector_d  scaled_vec = m_k.scaled_vector_d_object();

    Tangent_space_basis ts(i);
    ts.reserve(m_intrinsic_dim);
    ts.push_back(scaled_vec(t, FT(1)/CGAL::sqrt(sqlen(t))));
    if (i != std::numeric_limits<std::size_t>::max())
      m_are_tangent_spaces_computed[i] = true;

    return ts;

#elif defined(CGAL_MESH_D_COMPUTE_TANGENT_PLANES_FOR_SPHERE_3)

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

    if (i != std::numeric_limits<std::size_t>::max())
      m_are_tangent_spaces_computed[i] = true;

    return ts;

#elif defined(CGAL_MESH_D_COMPUTE_TANGENT_PLANES_FOR_TORUS_D)

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

    if (i != std::numeric_limits<std::size_t>::max())
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

#ifdef CGAL_MESH_D_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
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
#ifdef CGAL_MESH_D_ADD_NOISE_TO_TANGENT_SPACE
        mat_points(j, i) += m_random_generator.get_double(
            -0.5*m_half_sparsity, 0.5*m_half_sparsity);
#endif
      }
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    Tangent_space_basis tsb(i);

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

    if (i != std::numeric_limits<std::size_t>::max())
      m_are_tangent_spaces_computed[i] = true;

    //*************************************************************************

    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, FT(1)/sqrt(sqlen(n)));
    //std::cerr << "IP = " << scalar_pdct(n, ts[0]) << " & " << scalar_pdct(n, ts[1]) << std::endl;

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

  Vertex_set get_neighbor_vertices(
    Triangulation const& tr, Tr_vertex_handle vh, 
    bool keep_infinite_vertex = false)
  {
    Vertex_set neighbors;
    std::vector<Tr_full_cell_handle> incident_cells;
    tr.incident_full_cells(
      vh, std::back_inserter(incident_cells));

    typename std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
    typename std::vector<Tr_full_cell_handle>::const_iterator it_c_end = incident_cells.end();
    // For each cell
    for (; it_c != it_c_end ; ++it_c)
    {
      for (int j = 0 ; j < m_ambient_dim + 1 ; ++j)
      {
        Tr_vertex_handle v = (*it_c)->vertex(j);
        if (keep_infinite_vertex
          || v->data() != std::numeric_limits<std::size_t>::max())
        {
          neighbors.insert(v);
        }
      }
    }
    neighbors.erase(vh);
    return neighbors;
  }

  template <typename Vertex_range>
  Vertex_set get_common_neighbor_vertices(
    Triangulation const& tr, Vertex_range const& vertices)
  {
    Vertex_range::const_iterator vh_it = vertices.begin();
    Vertex_set common_neighbors =
      get_neighbor_vertices(tr, *vh_it);
    ++vh_it;

    for (; vh_it != vertices.end() ; ++vh_it)
    {
      Vertex_set neighbors = get_neighbor_vertices(tr, *vh_it);
      Vertex_set former_common_neighbors = common_neighbors;
      common_neighbors.clear();
      std::set_intersection(
        former_common_neighbors.begin(), former_common_neighbors.end(), 
        neighbors.begin(), neighbors.end(),
        std::inserter(common_neighbors, common_neighbors.begin()));
    }

    return common_neighbors;
  }

  template <typename Basis_, typename Projected_point, typename Projection_space_kernel>
  Point unproject_point(Projected_point const& lp,
    Basis_ const& basis, Projection_space_kernel const& lk,
    Point const* origin = NULL) const
  {
    typename K::Translated_point_d k_transl =
      m_k.translated_point_d_object();
    typename K::Scaled_vector_d k_scaled_vec =
      m_k.scaled_vector_d_object();
    typename Projection_space_kernel::Compute_coordinate_d coord =
      lk.compute_coordinate_d_object();

    Point global_point = (origin ? *origin : m_points[basis.origin()]);
    for (int i = 0 ; i < basis.dimension() ; ++i)
    {
      global_point = k_transl(
        global_point, k_scaled_vec(basis[i], coord(lp, i)));
    }

    return global_point;
  }

  // Project the point on a subspace
  // Resulting point coords are expressed in basis' space
  template <typename Projection_space_kernel, typename Basis_>
  typename Projection_space_kernel::Point_d project_point(
    const Point &p, const Basis_ &basis, Projection_space_kernel const& psk,
    Point const* origin = NULL) const
  {
    typename K::Scalar_product_d scalar_pdct =
      m_k.scalar_product_d_object();
    typename K::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();

    Point const& orig = (origin ? *origin : m_points[basis.origin()]);
    Vector v = diff_points(p, orig);

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    coords.reserve(basis.dimension());
    for (std::size_t i = 0 ; i < basis.dimension() ; ++i)
    {
      // Local coords are given by the scalar product with the vectors of basis
      FT coord = scalar_pdct(v, basis[i]);
      coords.push_back(coord);
    }

    return psk.construct_point_d_object()(static_cast<int>(
      coords.size()), coords.begin(), coords.end());
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  // Resulting point coords are expressed in tsb's space
  template <typename Projection_space_kernel, typename Basis_>
  typename Projection_space_kernel::Weighted_point_d 
  project_point_and_compute_weight(
    const Point &p, const FT w, const Basis_ &basis, 
    Projection_space_kernel const& psk,
    Point const* origin = NULL) const
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

    typename Projection_space_kernel::Construct_point_d local_constr_pt =
      psk.construct_point_d_object();

    Point const& orig = (origin ? *origin : m_points[basis.origin()]);
    Vector v = diff_points(p, orig);

    // Same dimension? Then the weight is 0
    bool same_dim = (point_dim == basis.dimension());

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(ccci(orig), ccci(orig, 0));
    coords.reserve(basis.dimension());
    for (std::size_t i = 0 ; i < basis.dimension() ; ++i)
    {
      // Local coords are given by the scalar product with the vectors of basis
      FT c = scalar_pdct(v, basis[i]);
      coords.push_back(c);

      // p_proj += c * basis[i]
      for (int j = 0 ; j < point_dim ; ++j)
        p_proj[j] += c * coord(basis[i], j);
    }
    
    // Same dimension? Then the weight is 0
    FT sq_dist_to_proj_pt = 0;
    if (!same_dim)
    {
      Point projected_pt = constr_pt(point_dim, p_proj.begin(), p_proj.end());
      sq_dist_to_proj_pt = m_k.squared_distance_d_object()(p, projected_pt);
    }

    return psk.construct_weighted_point_d_object()
    (
      local_constr_pt(
        static_cast<int>(coords.size()), coords.begin(), coords.end()),
      w - sq_dist_to_proj_pt
    );
  }

  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  // Note that the computation is made in global coordinates.
  template <typename Point_range_a, typename Point_range_b>
  CGAL::Quadratic_program_solution<ET>
    compute_voronoi_face_and_tangent_subspace_LP_problem(
    Point_range_a const& P,
    Point_range_b const& Q,
    Orthogonal_space_basis const& orthogonal_subspace_basis) const
  {
    // Notations:
    // Fv: Voronoi k-face
    // Fd: dual, (D-k)-face of Delaunay (p0, p1, ..., pn)

    typename K::Scalar_product_d scalar_pdct = m_k.scalar_product_d_object();
    typename K::Point_to_vector_d pt_to_vec = m_k.point_to_vector_d_object();
    typename K::Compute_coordinate_d coord = m_k.compute_coordinate_d_object();

    std::size_t card_P = P.size();
    std::size_t card_Q = Q.size();

    // Linear solver
    typedef CGAL::Quadratic_program<FT> Linear_program;
    typedef CGAL::Quadratic_program_solution<ET> LP_solution;

    Linear_program lp(CGAL::SMALLER, false);
    int current_row = 0;

    //=========== First set of equations ===========
    // For points pi in P
    //   2(p0 - pi).x = p0^2 - pi^2
    typename Point_range_a::const_iterator it_p = P.begin();
    Point const& p0 = *it_p;
    FT p0_dot_p0 = scalar_pdct(pt_to_vec(p0), pt_to_vec(p0));
    ++it_p;
    for (typename Point_range_a::const_iterator it_p_end = P.end() ;
      it_p != it_p_end ; ++it_p)
    {
      Point const& pi = *it_p;

      for (int k = 0 ; k < m_ambient_dim ; ++k)
        lp.set_a(k, current_row, 2 * (coord(p0, k) - coord(pi, k)));

      FT pi_dot_pi = scalar_pdct(pt_to_vec(pi), pt_to_vec(pi));
      lp.set_b(current_row, p0_dot_p0 - pi_dot_pi);
      lp.set_r(current_row, CGAL::EQUAL);

      ++current_row;
    }

    //=========== Second set of equations ===========
    // For each point qi in Q
    //  2(qi - p0).x <= qi^2 - p0^2
    for (typename Point_range_b::const_iterator it_q = Q.begin(),
      it_q_end = Q.end() ;
      it_q != it_q_end ; ++it_q)
    {
      Point const& qi = *it_q;

      for (int k = 0 ; k < m_ambient_dim ; ++k)
        lp.set_a(k, current_row, 2 * (coord(qi, k) - coord(p0, k)));

      FT qi_dot_qi = scalar_pdct(pt_to_vec(qi), pt_to_vec(qi));
      lp.set_b(current_row, qi_dot_qi - p0_dot_p0);

      ++current_row;
    }

    //=========== Third set of equations ===========
    // For each vector bi of OSB, (x-p).bi = 0
    //   <=> bi.x =  bi.p
    for (Orthogonal_space_basis::const_iterator it_osb =
      orthogonal_subspace_basis.begin(),
      it_osb_end = orthogonal_subspace_basis.end() ;
    it_osb != it_osb_end ; ++it_osb)
    {
      Vector const& bi = *it_osb;

      for (int k = 0 ; k < m_ambient_dim ; ++k)
      {
        lp.set_a(k, current_row, coord(bi, k));
      }

      FT bi_dot_p = scalar_pdct(bi,
        pt_to_vec(m_points[orthogonal_subspace_basis.origin()]));
      lp.set_b(current_row, bi_dot_p);
      lp.set_r(current_row, CGAL::EQUAL);

      ++current_row;
    }

    //=========== Other LP parameters ===========
    lp.set_c(0, 1); // Minimize x[0]

    //=========== Solve =========================
    LP_solution solution = CGAL::solve_linear_program(lp, ET());
    return solution;
  }

  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  template <typename Point_range_a, typename Point_range_b>
  bool does_voronoi_face_and_tangent_subspace_intersect(
    Point_range_a const& P,
    Point_range_b const& Q,
    Orthogonal_space_basis const& orthogonal_subspace_basis) const
  {
    return compute_voronoi_face_and_tangent_subspace_LP_problem(
      P, Q, orthogonal_subspace_basis).status() == CGAL::QP_OPTIMAL;
  }


  // Returns any point of the intersection between a Voronoi cell and a
  // tangent space.
  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Q: vertices which are common neighbors of all vertices of P
  // Return value: the point coordinates are expressed in global coordinates
  template <typename Point_range_a, typename Point_range_b>
  boost::optional<Point>
    compute_voronoi_face_and_tangent_subspace_intersection(
    Point_range_a const& P,
    Point_range_b const& Q,
    Orthogonal_space_basis const& orthogonal_subspace_basis) const
  {
    typedef CGAL::Quadratic_program_solution<ET> LP_solution;

    LP_solution sol = compute_voronoi_face_and_tangent_subspace_LP_problem(
      P, Q, orthogonal_subspace_basis);

    boost::optional<Point> ret;
    if (sol.status() == CGAL::QP_OPTIMAL)
    {
      std::vector<FT> p;
      p.reserve(m_ambient_dim);
      for (LP_solution::Variable_value_iterator
        it_v = sol.variable_values_begin(),
        it_v_end = sol.variable_values_end() ;
      it_v != it_v_end ; ++it_v)
      {
        p.push_back(to_double(*it_v));
      }
      CGAL_assertion(p.size() == m_ambient_dim);
      ret = m_k.construct_point_d_object()(m_ambient_dim, p.begin(), p.end());
    }
    else
    {
      ret = boost::none;
    }

    return ret;
  }

#ifdef CGAL_MESH_D_USE_LINEAR_PROG_TO_COMPUTE_INTERSECTION
  // CJTODO TEMP: this is the old slow code => remove this macro

  // Returns any point of the intersection between aff(voronoi_cell) and a
  // tangent space.
  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Return value: the point coordinates are expressed in the tsb base
  template <typename Point_range_a>
  Point
  compute_aff_of_voronoi_face_and_tangent_subspace_intersection(
    Point_range_a const& P,
    Orthogonal_space_basis const& orthogonal_subspace_basis) const
  {
    // As we're only interested by aff(v), Q is empty
    return *compute_voronoi_face_and_tangent_subspace_intersection(
      P, std::vector<typename Point_range_a::value_type>(),
      orthogonal_subspace_basis);
  }

#else 

  // Returns any point of the intersection between aff(voronoi_cell) and a
  // tangent space.
  // P: dual face in Delaunay triangulation (p0, p1, ..., pn)
  // Return value: the point coordinates are expressed in the tsb base
  template <typename Point_range>
  Point
    compute_aff_of_voronoi_face_and_tangent_subspace_intersection(
    Point_range const& P,
    Tangent_space_basis const& tangent_subspace_basis,
    Point const* origin = NULL) const
  {
    std::vector<Local_weighted_point> projected_pts;
    for (auto const& p : P)
    {
      projected_pts.push_back(project_point_and_compute_weight(
        p, FT(0), tangent_subspace_basis, m_lk, origin));
    }

    Local_weighted_point power_center =
      m_lk.power_center_d_object()(projected_pts.begin(), projected_pts.end());

    return unproject_point(
      m_lk.point_drop_weight_d_object()(power_center), 
      tangent_subspace_basis, m_lk, origin);
  }

#endif

  std::map<Vertex_set, Vertex_set> get_k_faces_and_neighbor_vertices(Triangulation const& tr)
  {
    typedef std::map<Vertex_set, Vertex_set> K_faces_and_neighbor_vertices;

    // Map that associate a k-face F and the points of its k+1-cofaces
    // (except the points of F). Those points are called its "neighbors".
    K_faces_and_neighbor_vertices faces_and_neighbors;

    // Fill faces_and_neighbors
    typedef K_faces_and_neighbor_vertices::const_iterator FaN_it;
    // Parse cells
    for (Tr_finite_full_cell_const_iterator cit = tr.finite_full_cells_begin() ;
      cit != tr.finite_full_cells_end() ; ++cit)
    {
      // Add each k-face to faces_and_neighbors
      std::vector<bool> booleans(m_ambient_dim + 1, true);
      std::fill(
        booleans.begin(),
        booleans.begin() + m_ambient_dim - m_intrinsic_dim,
        false);
      do
      {
        Vertex_set k_face;
        std::vector<Tr_vertex_handle> remaining_vertices;
        for (int i = 0 ; i < m_ambient_dim + 1 ; ++i)
        {
          if (booleans[i])
            k_face.insert(cit->vertex(i));
          else
            remaining_vertices.push_back(cit->vertex(i));
        }

        faces_and_neighbors[k_face].insert(
          remaining_vertices.begin(), remaining_vertices.end());

      } while (std::next_permutation(booleans.begin(), booleans.end()));
    }

    return faces_and_neighbors;
  }

  // Returns the dimension of the ith local triangulation
  // This is particularly useful for the alpha-TC
  int tangent_basis_dim(std::size_t i) const
  {
    return m_tangent_spaces[i].dimension();
  }

private:
  std::ostream &export_vertices_to_off(
    std::ostream & os, std::size_t &num_vertices) const
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

    int num_coords = min(m_ambient_dim, 3);
#ifdef CGAL_MESH_D_EXPORT_NORMALS
    OS_container::const_iterator it_os = m_orth_spaces.begin();
#endif
    typename Points::const_iterator it_p = m_points.begin();
    typename Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for (std::size_t i = 0 ; it_p != it_p_end ; ++it_p, ++i)
    {
      Point const& p = *it_p;
      for (int ii = 0 ; ii < N ; ++ii)
      {
        int i = 0;
#if BETTER_EXPORT_FOR_FLAT_TORUS
        // For flat torus
        os << (2 + 1 * CGAL::to_double(coord(p, 0))) * CGAL::to_double(coord(p, 2)) << " "
           << (2 + 1 * CGAL::to_double(coord(p, 0))) * CGAL::to_double(coord(p, 3)) << " "
           << 1 * CGAL::to_double(coord(p, 1));
#else
        for ( ; i < num_coords ; ++i)
          os << CGAL::to_double(coord(p, i)) << " ";
#endif
        if (i == 2)
          os << "0";

#ifdef CGAL_MESH_D_EXPORT_NORMALS
        for (i = 0 ; i < num_coords ; ++i)
          os << " " << CGAL::to_double(coord(*it_os->begin(), i));
#endif
        os << std::endl;
      }
#ifdef CGAL_MESH_D_EXPORT_NORMALS
      ++it_os;
#endif
    }

    num_vertices = N*m_points.size();
    return os;
  }

  template <typename Indexed_simplex_range = void>
  std::ostream &export_simplices_to_off(
    std::ostream & os, std::size_t &num_OFF_simplices,
    Indexed_simplex_range const *p_simpl_to_color_in_red = NULL,
    Indexed_simplex_range const *p_simpl_to_color_in_green = NULL,
    Indexed_simplex_range const *p_simpl_to_color_in_blue = NULL)
    const
  {
    typedef Simplicial_complex::Simplex                     Simplex;
    typedef Simplicial_complex::Simplex_set                 Simplex_set;

    // If m_intrinsic_dim = 1, each point is output two times
    // (see export_vertices_to_off)
    num_OFF_simplices = 0;
    std::size_t num_maximal_simplices = 0;

    typename Simplex_set::const_iterator it_s =
      m_complex.simplex_range().begin();
    typename Simplex_set::const_iterator it_s_end =
      m_complex.simplex_range().end();
    // For each simplex
    for (; it_s != it_s_end ; ++it_s)
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
        for (; it != c.end() ; ++it)
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
        } while (std::next_permutation(booleans.begin(), booleans.end()));
      }

      // For each cell
      Triangles::const_iterator it_tri = triangles.begin();
      Triangles::const_iterator it_tri_end = triangles.end();
      for (; it_tri != it_tri_end ; ++it_tri)
      {
        // Don't export infinite cells
        if (is_infinite(*it_tri))
          continue;

        os << 3 << " ";
        Indexed_simplex::const_iterator it_point_idx = it_tri->begin();
        for (; it_point_idx != it_tri->end() ; ++it_point_idx)
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
        os << std::endl;
      }
    }

#ifdef CGAL_MESH_D_VERBOSE
    std::cerr << std::endl
      << "=========================================================="
      << std::endl
      << "Export from complex to OFF:\n"
      << "  * Number of vertices: " << m_points.size() << std::endl
      << "  * Total number of maximal simplices: " << num_maximal_simplices
      << std::endl
      << "=========================================================="
      << std::endl;
#endif

    return os;
  }

public: 
  template <typename Indexed_simplex_range = void>
  std::ostream &export_to_off(
    std::ostream & os,
    Indexed_simplex_range const *p_simpl_to_color_in_red = NULL,
    Indexed_simplex_range const *p_simpl_to_color_in_green = NULL,
    Indexed_simplex_range const *p_simpl_to_color_in_blue = NULL) const
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

    if (m_intrinsic_dim < 1 || m_intrinsic_dim > 3)
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
      output, num_simplices, p_simpl_to_color_in_red,
      p_simpl_to_color_in_green, p_simpl_to_color_in_blue);

#ifdef CGAL_MESH_D_EXPORT_NORMALS
    os << "N";
#endif

    os << "OFF \n"
      << num_vertices << " "
      << num_simplices << " "
      << "0 \n"
      << output.str();

    return os;
  }

private:
  const K                   m_k;
  const LK                  m_lk;
  const int                 m_intrinsic_dim;
  const double              m_half_sparsity;
  const double              m_sq_half_sparsity;
  const int                 m_ambient_dim;

  Points                    m_points;

  Points_ds                 m_points_ds;
  std::vector<bool>         m_are_tangent_spaces_computed;
  TS_container              m_tangent_spaces;
  OS_container              m_orth_spaces;

#ifdef CGAL_MESH_D_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  Points                    m_points_for_tse;
  Points_ds                 m_points_ds_for_tse;
#endif

  mutable CGAL::Random      m_random_generator;

  Simplicial_complex        m_complex;

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // MESH_D_H
