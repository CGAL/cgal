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

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Tangential_complex/utilities.h>

#include <CGAL/Mesh_3/Profiling_tools.h>

#include <CGAL/IO/Triangulation_off_ostream.h> // CJTODO TEMP

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <vector>
#include <utility>
#include <sstream>
#include <iostream>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/parallel_for.h>
#endif

namespace CGAL {

using namespace Tangential_complex_;
  
/// The class Tangential_complex represents a tangential complex
template <
  typename Kernel, 
  int Intrinsic_dimension,
  typename Concurrency_tag = CGAL::Parallel_tag,
  typename Tr = Regular_triangulation
  <
    Regular_triangulation_euclidean_traits<
      Epick_d<Dimension_tag<Intrinsic_dimension> > >,

    Triangulation_data_structure
    <
      typename Regular_triangulation_euclidean_traits<
        Epick_d<Dimension_tag<Intrinsic_dimension> > >::Dimension,
      Triangulation_vertex<Regular_triangulation_euclidean_traits<
        Epick_d<Dimension_tag<Intrinsic_dimension> > >, std::size_t >,
      Triangulation_full_cell<Regular_triangulation_euclidean_traits<
        Epick_d<Dimension_tag<Intrinsic_dimension> > > >
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
  typedef typename Triangulation::Point               Tr_point;
  typedef typename Triangulation::Bare_point          Tr_bare_point;
  typedef typename Triangulation::Vertex_handle       Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle    Tr_full_cell_handle;
  
  typedef typename std::vector<Vector>                Tangent_space_basis;

  typedef std::pair<Triangulation*, Tr_vertex_handle> Tr_and_VH;
  typedef typename std::vector<Point>                 Point_container;
  typedef typename std::vector<Tr_and_VH>             Tr_container;
  typedef typename std::vector<Tangent_space_basis>    TS_container;

  // Stores the index of the original Point in the ambient space
  /*struct Tr_point_with_index 
  : public Tr_point
  {
    Tr_point_with_index(const Tr_point &p, std::size_t i)
      : Tr_point(p), index(i) {}

    std::size_t index;
  };*/

public:
  /// Constructor
  Tangential_complex(const Kernel &k = Kernel())
  : m_k(k){}
  
  /// Constructor for a range of points
  template <typename InputIterator>
  Tangential_complex(InputIterator first, InputIterator last, 
                     const Kernel &k = Kernel())
  : m_k(k), m_points(first, last) {}

  /// Destructor
  ~Tangential_complex() {}

  void compute_tangential_complex()
  {
#ifdef CGAL_TC_PROFILING
    WallClockTimer t;
#endif

    // We need to do that because we don't want the container to copy the
    // already-computed triangulations (while resizing) since it would
    // invalidate the vertex handles stored beside the triangulations
    m_triangulations.resize(
      m_points.size(), 
      std::make_pair((Triangulation*)NULL, Tr_vertex_handle()));
    m_tangent_spaces.resize(m_points.size());
    
#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      // Apply moves in triangulation
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

  std::ostream &export_to_off(std::ostream & os)
  {
    const int ambient_dim = Ambient_dimension<Point>::value;
    if (ambient_dim < 2 || ambient_dim > 3)
    {
      std::cerr << "Error: export_to_off => ambient dimension should be 2 or 3.";
      os << "Error: export_to_off => ambient dimension should be 2 or 3.";
      return os;
    }

    if (Intrinsic_dimension < 1 || Intrinsic_dimension > 3)
    {
      std::cerr << "Error: export_to_off => intrinsic dimension should be between 1 and 3.";
      os << "Error: export_to_off => intrinsic dimension should be between 1 and 3.";
      return os;
    }

    std::stringstream output;

    //******** VERTICES ************

    Point_container::const_iterator it_p = m_points.begin();
    Point_container::const_iterator it_p_end = m_points.end();
    // For each point p
    for ( ; it_p != it_p_end ; ++it_p)
    {
      int i = 0;
      for ( ; i < ambient_dim ; ++i)
        output << (*it_p)[i] << " ";
      if (i == 2)
        output << "0";
      output << std::endl;
    }

    //******** CELLS ************

    std::size_t num_cells = 0;
    Tr_container::const_iterator it_tr = m_triangulations.begin();
    Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for ( ; it_tr != it_tr_end ; ++it_tr)
    {
      const Triangulation &tr = *it_tr->first;
      Tr_vertex_handle center_vh = it_tr->second;

      std::vector<Tr_full_cell_handle> incident_cells;
      tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

      std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
      std::vector<Tr_full_cell_handle>::const_iterator it_c_end= incident_cells.end();
      // For each triangulation
      for ( ; it_c != it_c_end ; ++it_c)
      {
        output << Intrinsic_dimension + 1 << " ";
        for (int i = 0 ; i < Intrinsic_dimension + 1 ; ++i)
          output << (*it_c)->vertex(i)->data() << " ";
        output << std::endl;
        ++num_cells;
      }
    }
    
    os << "OFF \n"
       << m_points.size() << " " 
       << num_cells << " "
       << "0 \n"
       << output.str();

    return os;
  }

private:

  class Compare_distance_to_ref_point
  {
  public:
    Compare_distance_to_ref_point(Point const& ref, Kernel const& k)
      : m_ref(ref), m_k(k) {}

    bool operator()(Point const& p1, Point const& p2)
    {
      Kernel::Squared_distance_d sqdist = m_k.squared_distance_d_object();
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
    Compute_tangent_triangulation(Tangential_complex &tc)
    : m_tc(tc)
    {}

    // Constructor
    Compute_tangent_triangulation(const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
        m_tc.compute_tangent_triangulation(i);
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  void compute_tangent_triangulation(std::size_t i)
  {
    Triangulation *p_local_tr =
      m_triangulations[i].first = 
        new Triangulation(Intrinsic_dimension);
    const Tr_traits &local_tr_traits = p_local_tr->geom_traits();
    Tr_vertex_handle &center_vertex = m_triangulations[i].second;

    // Estimate the tangent space
    const Point &center_pt = m_points[i];
    m_tangent_spaces[i] = compute_tangent_space(center_pt);
      
    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p)
    //***************************************************

    // First, compute the projected points
    std::vector<Tr_point> projected_points;
    FT max_squared_weight = 0;
    projected_points.reserve(m_points.size() - 1);
    Point_container::const_iterator it_p = m_points.begin();
    Point_container::const_iterator it_p_end = m_points.end();
    for (std::size_t j = 0 ; it_p != it_p_end ; ++it_p, ++j)
    {
      // ith point = p, which is already inserted
      if (j != i)
      {
        Tr_point wp = project_point(*it_p, center_pt, m_tangent_spaces[i]);
        projected_points.push_back(wp);
        FT w = local_tr_traits.point_weight_d_object()(wp);
        if (w > max_squared_weight)
          max_squared_weight = w;
      }
    }

    // Now we can insert the points

    // Insert p
    Tr_point wp = local_tr_traits.construct_weighted_point_d_object()(
      local_tr_traits.construct_point_d_object()(0, 0),
      CGAL::sqrt(max_squared_weight));
    center_vertex = p_local_tr->insert(wp);
    center_vertex->data() = i;
    //std::cerr << "Inserted CENTER POINT of weight " << CGAL::sqrt(max_squared_weight) << std::endl;
      
    /*std::cerr << 0 << " "
              << 0 << " "
              << CGAL::sqrt(max_squared_weight) << std::endl;*/

    // Insert the other points
    std::vector<Tr_point>::const_iterator it_wp = projected_points.begin();
    it_p = m_points.begin();
    for (std::size_t j = 0 ; it_p != it_p_end ; ++it_p, ++j)
    {
      // ith point = p, which is already inserted
      if (j != i)
      {
        // CJTODO TEMP: for test only
        /*if (local_tr_traits.squared_distance_d_object()(
          local_tr_traits.point_drop_weight_d_object()(wp), 
          local_tr_traits.point_drop_weight_d_object()(*it_wp)) > 1)
        {
          ++it_wp;
          continue;
        }*/

        FT squared_dist_to_tangent_plane = 
          local_tr_traits.point_weight_d_object()(*it_wp);
        FT w = CGAL::sqrt(max_squared_weight - squared_dist_to_tangent_plane);
        Tr_point wp = local_tr_traits.construct_weighted_point_d_object()(
          local_tr_traits.point_drop_weight_d_object()(*it_wp),
          w);
        /*Tr_bare_point bp = traits.point_drop_weight_d_object()(*it_wp);
        Tr_point wp(traits.point_drop_weight_d_object()(*it_wp), w);*/
          
        Tr_vertex_handle vh = p_local_tr->insert_if_in_star(wp, center_vertex);
        //Tr_vertex_handle vh = p_local_tr->insert(wp);
        if (vh != Tr_vertex_handle())
        {
          /*std::cerr << traits.point_drop_weight_d_object()(*it_wp)[0] << " "
                    << traits.point_drop_weight_d_object()(*it_wp)[1] << " "
                    << w << std::endl;*/
          vh->data() = j;
        }
        ++it_wp;
      }
    }

    // CJTODO DEBUG
    //std::cerr << "\nChecking topology and geometry..."
    //          << (p_local_tr->is_valid(true) ? "OK.\n" : "Error.\n");
    // DEBUG: output the local mesh into an OFF file
    //std::stringstream sstr;
    //sstr << "data/local_tri_" << i << ".off";
    //std::ofstream off_stream_tr(sstr.str());
    //CGAL::export_triangulation_to_off(off_stream_tr, *p_local_tr);
  }

  Tangent_space_basis compute_tangent_space(const Point &p) const
  {
    // Kernel functors
    Kernel::Construct_vector_d      constr_vec = m_k.construct_vector_d_object();
    Kernel::Squared_length_d        sqlen      = m_k.squared_length_d_object();
    Kernel::Scaled_vector_d         scaled_vec = m_k.scaled_vector_d_object();
    //Kernel::Scalar_product_d        inner_pdct = m_k.scalar_product_d_object();
    //Kernel::Difference_of_vectors_d diff_vec   = m_k.difference_of_vectors_d_object();
    Get_functor<Kernel, Scalar_product_tag>::type inner_pdct(m_k); // CJTODO TEMP
    Get_functor<Kernel, Difference_of_vectors_tag>::type diff_vec(m_k);

    // CJTODO: do better than that (ANN?)
    typedef std::set<Point, Compare_distance_to_ref_point> Sorted_points;
    Sorted_points sorted_points(
      Compare_distance_to_ref_point(p, m_k));
    sorted_points.insert(m_points.begin(), m_points.end());

    //******************************* PCA *************************************

    const int amb_dim = Ambient_dimension<Point>::value;
    Eigen::MatrixXd mat(NUM_POINTS_FOR_PCA, amb_dim);
    int j = 0;
    for (Sorted_points::const_iterator it = sorted_points.begin() ; 
         j < NUM_POINTS_FOR_PCA ; ++it, ++j)
    {
      for (int i = 0 ; i < amb_dim ; ++i)
        mat(j, i) = (*it)[i];
    }
    Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    Tangent_space_basis ts;
    for (int i = amb_dim - 1 ; i >= amb_dim - Intrinsic_dimension ; --i)
    {
      ts.push_back(constr_vec(
        amb_dim, 
        eig.eigenvectors().col(i).data(), 
        eig.eigenvectors().col(i).data() + amb_dim));
    }

    //*************************************************************************

    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, 1./sqrt(sqlen(n)));
    //std::cerr << "IP = " << inner_pdct(n, ts[0]) << " & " << inner_pdct(n, ts[1]) << std::endl;

    return compute_gram_schmidt_basis(ts, m_k);

    /*
    // CJTODO: this is only for a sphere in R^3
    Vector t1(-p[1] - p[2], p[0], p[0]);
    Vector t2(p[1] * t1[2] - p[2] * t1[1],
              p[2] * t1[0] - p[0] * t1[2],
              p[0] * t1[1] - p[1] * t1[0]);
    
    // Normalize t1 and t2
    Get_functor<Kernel, Scaled_vector_tag>::type scale(m_k);

    Tangent_space_basis ts;
    ts.reserve(Intrinsic_dimension);
    ts.push_back(scale(t1, 1./CGAL::sqrt(sqlen(t1))));
    ts.push_back(scale(t2, 1./CGAL::sqrt(sqlen(t2))));

    return ts;

    // Alternative code (to be used later)
    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, 1./sqrt(sqlen(n)));
    //Vector t1(12., 15., 65.);
    //Vector t2(32., 5., 85.);
    //Tangent_space_basis ts;
    //ts.reserve(Intrinsic_dimension);
    //ts.push_back(diff_vec(t1, scaled_vec(n, inner_pdct(t1, n))));
    //ts.push_back(diff_vec(t2, scaled_vec(n, inner_pdct(t2, n))));
    //return compute_gram_schmidt_basis(ts, m_k);
    */
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_point project_point(const Point &p, const Point &origin, 
                         const Tangent_space_basis &ts) const
  {
    Get_functor<Kernel, Scalar_product_tag>::type inner_pdct(m_k);
    Get_functor<Kernel, Difference_of_points_tag>::type diff_points(m_k);
  
    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(origin.cartesian_begin(), origin.cartesian_end());
    coords.reserve(Intrinsic_dimension);
    for (std::size_t i = 0 ; i < Intrinsic_dimension ; ++i)
    {
      // Compute the inner product p * ts[i]
      Vector v = diff_points(p, origin);
      FT coord = inner_pdct(v, ts[i]);
      coords.push_back(coord);

      // p_proj += coord * v;
      for (int j = 0 ; j < Ambient_dimension<Point>::value ; ++j)
        p_proj[i] += coord * ts[i][j];
    }

    Point projected_pt(Ambient_dimension<Point>::value, 
                       p_proj.begin(), p_proj.end());
    return Tr_point(
      Tr_bare_point(Intrinsic_dimension, coords.begin(), coords.end()), 
      m_k.squared_distance_d_object()(p, projected_pt));
  }

private:
  const Kernel        m_k;
  Point_container     m_points;
  TS_container        m_tangent_spaces;
  Tr_container        m_triangulations; // Contains the triangulations 
                                        // and their center vertex

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
